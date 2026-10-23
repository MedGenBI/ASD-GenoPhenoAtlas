import argparse
import math
import os
import re
from pathlib import Path
from statistics import NormalDist

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter


SCRIPT_NAME = "08_patient_level_eval"
TOP_NS = (1, 5, 10, 20)
METHOD_ORDER = ("GPS", "CADA", "Phen2Gene", "PhenoGenius", "GADO")
METHOD_COLORS = {
    "GPS": "#E4BFCD",
    "CADA": "#C6CCDD",
    "Phen2Gene": "#A8B2CB",
    "PhenoGenius": "#818EAD",
    "GADO": "#657493",
}
METHOD_EDGES = {
    "GPS": "#8F5A6D",
    "CADA": "#ACB4C5",
    "Phen2Gene": "#98A1BA",
    "PhenoGenius": "#666F87",
    "GADO": "#4E5C78",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate patient-level Top-N accuracy and optional tool comparison."
    )
    parser.add_argument(
        "--benchmark-summary",
        default=None,
        help=(
            "Optional sample-level benchmark Top-N summary TSV. GPS rows are ignored; "
            "non-GPS methods are compared with the supplied GPS prediction."
        ),
    )
    parser.add_argument(
        "--gps-prediction",
        required=True,
        help="06/07 patient-level GPS prediction TSV used to calculate GPS results.",
    )
    parser.add_argument(
        "--analysis-label",
        required=True,
        help="Output label, for example derivation or external.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Project output root; results are written below 08_patient_level_eval.",
    )
    return parser.parse_args()


def canonical_method(value):
    text = str(value).strip()
    aliases = {
        "gps": "GPS",
        "cada": "CADA",
        "phen2gene": "Phen2Gene",
        "phenogenius": "PhenoGenius",
        "gado": "GADO",
    }
    return aliases.get(text.lower(), text)


def parse_top_n(value):
    if pd.isna(value):
        return None
    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return None
    numeric = float(numeric)
    return int(numeric) if numeric.is_integer() else None


def wilson_interval(hits, n, confidence=0.95):
    if n <= 0:
        return math.nan, math.nan
    z = NormalDist().inv_cdf(1.0 - (1.0 - confidence) / 2.0)
    p = hits / n
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    radius = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))
    radius /= denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def load_summary(path):
    frame = pd.read_csv(path, sep="\t")
    if frame.empty:
        raise ValueError(f"Comparison summary is empty: {path}")

    required = {"tool", "top_n", "n_evaluated", "topn_hits"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            "The comparison summary is missing required columns: "
            + ", ".join(missing)
        )

    records = []
    for _, row in frame.iterrows():
        method = canonical_method(row["tool"])
        if method == "GPS":
            continue

        top_n = parse_top_n(row["top_n"])
        if top_n not in TOP_NS:
            continue

        n = pd.to_numeric(row["n_evaluated"], errors="coerce")
        hits = pd.to_numeric(row["topn_hits"], errors="coerce")
        if (
            pd.isna(n)
            or pd.isna(hits)
            or n <= 0
            or hits < 0
            or hits > n
            or not float(n).is_integer()
            or not float(hits).is_integer()
        ):
            continue
        n = int(n)
        hits = int(hits)
        accuracy = hits / n

        lower, upper = wilson_interval(hits, n)
        records.append(
            {
                "method": method,
                "top_n": top_n,
                "hits": hits,
                "n": n,
                "accuracy": float(accuracy),
                "ci_lower": lower,
                "ci_upper": upper,
                "source": str(Path(path).resolve()),
            }
        )

    result = pd.DataFrame(records)
    if result.empty:
        raise ValueError(f"No sample-level Top-N rows for {TOP_NS}: {path}")
    return result


def parse_bool(value):
    text = str(value).strip().lower()
    return text in {"true", "1", "yes", "y", "t"}


def load_gps_prediction(path):
    frame = pd.read_csv(path, sep="\t")
    required = {"is_panel_carrier_sample", "best_true_rank"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            "GPS prediction TSV is missing required columns: " + ", ".join(missing)
        )

    frame = frame[frame["is_panel_carrier_sample"].map(parse_bool)].copy()
    ranks = pd.to_numeric(frame["best_true_rank"], errors="coerce")
    frame = frame.loc[ranks.notna()].copy()
    ranks = ranks.loc[frame.index]

    records = []
    for top_n in TOP_NS:
        hits = int((ranks <= top_n).sum())
        n = int(len(ranks))
        lower, upper = wilson_interval(hits, n)
        records.append(
            {
                "method": "GPS",
                "top_n": top_n,
                "hits": hits,
                "n": n,
                "accuracy": hits / n if n else math.nan,
                "ci_lower": lower,
                "ci_upper": upper,
                "source": str(Path(path).resolve()),
            }
        )
    return pd.DataFrame(records)


def combine_results(summary, gps_prediction):
    if gps_prediction is not None:
        summary = summary[summary["method"] != "GPS"].copy()
        summary = pd.concat([summary, gps_prediction], ignore_index=True)

    summary = summary[summary["method"].isin(METHOD_ORDER)].copy()
    summary = summary[summary["top_n"].isin(TOP_NS)].copy()
    summary["method"] = pd.Categorical(
        summary["method"], categories=METHOD_ORDER, ordered=True
    )
    summary = summary.sort_values(["top_n", "method"]).reset_index(drop=True)
    duplicate = summary.duplicated(["method", "top_n"], keep=False)
    if duplicate.any():
        raise ValueError(
            "Comparison summary must contain one sample-level row per method and Top-N."
        )
    return summary


def configure_plot_style():
    """Use the same typography and export settings as the original workflow."""
    requested_font = "Arial"
    try:
        findfont(FontProperties(family=requested_font), fallback_to_default=False)
        resolved_font = requested_font
    except (OSError, ValueError):
        resolved_font = "Liberation Sans"
        try:
            findfont(FontProperties(family=resolved_font), fallback_to_default=False)
        except (OSError, ValueError):
            resolved_font = "DejaVu Sans"

    plt.rcParams.update(
        {
            "font.family": resolved_font,
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Liberation Sans",
                "DejaVu Sans",
            ],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "text.color": "#222222",
            "axes.labelcolor": "#222222",
            "axes.edgecolor": "#333333",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def configure_facet_y_ticks(axis, y_min, y_max):
    """Match the original independent y-axis tick locator and formatter."""
    y_span = max(1e-12, float(y_max - y_min))
    axis.yaxis.set_major_locator(
        MaxNLocator(
            nbins=5,
            steps=[1, 2, 2.5, 5, 10],
            min_n_ticks=4,
            prune=None,
        )
    )
    approximate_step = y_span / 5.0
    decimals = 2 if approximate_step < 0.1 else 1
    axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals}f"))
    axis.tick_params(axis="y", labelsize=6.5)


def resolve_independent_facet_limits(panel):
    """Set the original data-driven limits for one Top-N facet."""
    data_min = float(panel["accuracy"].min())
    data_max = float(panel["accuracy"].max())
    padding = max(0.025, (data_max - data_min) * 0.18)
    y_min = max(0.0, math.floor((data_min - padding) / 0.05) * 0.05)
    y_max = min(1.0, math.ceil((data_max + padding) / 0.05) * 0.05)
    if y_max <= y_min:
        y_min = max(0.0, y_max - 0.1)
    return y_min, y_max


def style_reference_bar_axis(axis):
    """Apply the original boxed axes and tick styling."""
    axis.grid(False)
    for side in ("left", "bottom", "top", "right"):
        spine = axis.spines[side]
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(0.72)
    axis.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=2.6,
        width=0.7,
        color="#333333",
        pad=2.0,
    )


def annotate_accuracy(axis, bar, value, y_min, y_max, highlight):
    if not np.isfinite(value):
        return
    axis_span = max(0.05, y_max - y_min)
    offset = max(0.004, axis_span * 0.025)
    outside_y = value + offset
    near_upper_edge = outside_y > y_max - max(0.012, axis_span * 0.10)
    y_value = y_max - offset if near_upper_edge else outside_y
    axis.text(
        bar.get_x() + bar.get_width() / 2.0,
        y_value,
        f"{value:.3f}",
        ha="center",
        va="top" if near_upper_edge else "bottom",
        fontsize=6.0,
        fontweight="bold" if highlight else "normal",
        clip_on=True,
        zorder=5,
    )


def plot_faceted_bar(frame, output_dir, label):
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_style()
    methods = [method for method in METHOD_ORDER if method in set(frame["method"])]
    fig, axes = plt.subplots(1, len(TOP_NS), figsize=(7.09, 2.05), squeeze=False)
    axes = axes[0]

    for axis, top_n in zip(axes, TOP_NS):
        panel = frame[frame["top_n"] == top_n].set_index("method").reindex(methods)
        values = panel["accuracy"].to_numpy(dtype=float)
        positions = np.arange(len(methods))
        y_min, y_max = resolve_independent_facet_limits(panel)
        bars = axis.bar(
            positions,
            values,
            width=0.72,
            color=[METHOD_COLORS[m] for m in methods],
            edgecolor=[METHOD_EDGES[m] for m in methods],
            linewidth=0.65,
            zorder=2,
        )
        for method_index, (bar, value) in enumerate(zip(bars, values)):
            annotate_accuracy(
                axis,
                bar,
                value,
                y_min,
                y_max,
                highlight=method_index == 0,
            )
        axis.set_ylim(y_min, y_max)
        configure_facet_y_ticks(axis, y_min, y_max)
        axis.set_title(f"n={top_n}", fontsize=7, fontweight="bold", pad=3.5)
        axis.set_xticks(positions)
        axis.set_xticklabels(methods, rotation=43, ha="right")
        if axis is axes[0]:
            axis.set_ylabel("Top-n accuracy", fontsize=7)
        style_reference_bar_axis(axis)

    fig.subplots_adjust(
        left=0.075,
        right=0.995,
        top=0.84,
        bottom=0.31,
        wspace=0.34,
    )
    paths = []
    for suffix in ("pdf", "svg", "png"):
        path = output_dir / f"fig_topn_accuracy.{label}.{suffix}"
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if suffix == "png":
            kwargs["dpi"] = 600
        fig.savefig(path, **kwargs)
        paths.append(str(path))
    plt.close(fig)
    return paths


def main():
    args = parse_args()
    summary_path = None
    if args.benchmark_summary:
        summary_path = Path(args.benchmark_summary).expanduser().resolve()
        if not summary_path.is_file():
            raise FileNotFoundError(f"Comparison summary not found: {summary_path}")

    prediction_path = Path(args.gps_prediction).expanduser().resolve()
    if not prediction_path.is_file():
        raise FileNotFoundError(f"GPS prediction not found: {prediction_path}")

    analysis_dir = Path(args.outdir).expanduser().resolve() / SCRIPT_NAME
    output_dir = analysis_dir / "output"
    figure_dir = analysis_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading GPS patient-level prediction...")
    gps_prediction = load_gps_prediction(prediction_path)
    if summary_path is None:
        print("No benchmark summary supplied; writing GPS-only results.")
        summary = combine_results(gps_prediction, None)
    else:
        print("Loading non-GPS benchmark results...")
        summary = combine_results(load_summary(summary_path), gps_prediction)

    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", args.analysis_label.strip()) or "analysis"
    summary_path_out = output_dir / f"patient_level_topn_summary.{label}.tsv"
    summary.to_csv(summary_path_out, sep="\t", index=False)
    plot_faceted_bar(summary, figure_dir, label)

    print(f"Top-N summary written to: {summary_path_out}")
    print(f"Figure written to: {figure_dir}")


if __name__ == "__main__":
    main()
