import argparse
import math
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2


CODE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_ROOT not in sys.path:
    sys.path.insert(0, CODE_ROOT)

from src.config import OUTPUT_DIR


DEFAULT_TOP_PERCENTS = "0.2,0.1,0.05,0.01"
OUTPUT_FORMATS = ("png", "pdf", "svg")
PLOT_DPI = 600
FIGURE_WIDTH = 3.54
FIGURE_HEIGHT = 2.45
GRID_COLOR = "#E6EAF0"
TREND_LINE_COLOR = "#D2D2D2"
TREND_OR_COLOR = "#D66B72"


def read_table(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    separator = "," if str(path).lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=separator)


def read_gene_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Candidate-gene file not found: {path}")
    genes = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip().split("\t")[0].split(",")[0].strip()
            if gene and not gene.startswith("#"):
                genes.append(gene)
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("Candidate-gene file is empty.")
    return genes


def normalize_gene_column(frame, source_name):
    result = frame.copy()
    if "Gene" not in result.columns:
        aliases = {str(column).strip().lower(): column for column in result.columns}
        source = aliases.get("gene") or aliases.get("gene_symbol")
        if source is not None:
            result = result.rename(columns={source: "Gene"})
    if "Gene" not in result.columns:
        raise ValueError(f"{source_name} must contain a Gene column.")
    result["Gene"] = result["Gene"].astype(str).str.strip()
    return result[result["Gene"].ne("")].copy()


def fraction_value(value):
    text = str(value).strip()
    if text.endswith("%"):
        return float(text[:-1]) / 100.0
    number = float(text)
    return number / 100.0 if number > 1 else number


def fraction_label(value):
    return f"{100.0 * float(value):g}%"


def percentile_label(value):
    percentile = 100.0 * (1.0 - float(value))
    if np.isclose(percentile, round(percentile)):
        number = int(round(percentile))
        if 10 <= number % 100 <= 20:
            suffix = "th"
        else:
            suffix = {1: "st", 2: "nd", 3: "rd"}.get(number % 10, "th")
        return f"{number}{suffix}"
    return f"{percentile:g}th"


def percentile_tick_label(label):
    text = str(label).strip()
    for suffix in ("st", "nd", "rd", "th"):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    return rf"$>{text}^{{\mathrm{{th}}}}$"


def parse_top_percents(value):
    """Parse and validate the requested trend fractions."""
    fractions = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        fraction = fraction_value(token)
        if not np.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError(
                "--top-percents values must be finite and in the interval (0, 1]."
            )
        if fraction not in fractions:
            fractions.append(fraction)
    if not fractions:
        raise ValueError("--top-percents must contain at least one value.")
    return fractions


def or_statistics(row):
    """Compute uncorrected OR statistics from an integer 2 x 2 table."""
    count_names = ("n_top_carriers", "n_top_noncarriers",
                   "n_remaining_carriers", "n_remaining_noncarriers")
    values = []
    for name in count_names:
        try:
            value = float(pd.to_numeric(row.get(name), errors="coerce"))
        except (TypeError, ValueError):
            value = np.nan
        if not np.isfinite(value) or value < 0 or not value.is_integer():
            return {
                "or": np.nan,
                "or_ci_low": np.nan,
                "or_ci_upp": np.nan,
                "p": np.nan,
                "or_status": "invalid_counts",
            }
        values.append(int(value))

    tp, fp, fn, tn = values
    numerator = tp * tn
    denominator = fp * fn
    if numerator == 0:
        odds_ratio = 0.0
        status = "undefined_zero_over_zero" if denominator == 0 else "zero"
    elif denominator == 0:
        odds_ratio = np.inf
        status = "infinite"
    else:
        odds_ratio = float(numerator / denominator)
        status = "finite_positive"

    try:
        fisher_p = float(
            fisher_exact([[tp, fp], [fn, tn]], alternative="two-sided")[1]
        )
    except (ValueError, TypeError):
        fisher_p = np.nan

    if all(value > 0 for value in values) and np.isfinite(odds_ratio):
        ci_low, ci_upp = Table2x2(
            [[tp, fp], [fn, tn]], shift_zeros=False
        ).oddsratio_confint()
        ci_low = float(ci_low)
        ci_upp = float(ci_upp)
    else:
        ci_low = np.nan
        ci_upp = np.nan

    return {
        "or": odds_ratio,
        "or_ci_low": ci_low,
        "or_ci_upp": ci_upp,
        "p": fisher_p,
        "or_status": status,
    }


def select_candidate_panel(candidate_spec, available_genes, gene_metrics_file):
    """Select the fixed panel used by downstream patient-level analysis."""
    available = list(dict.fromkeys(str(gene) for gene in available_genes))
    available_set = set(available)
    spec = str(candidate_spec or "auto").strip() or "auto"

    if spec.lower() == "all":
        requested = sorted(available_set)
    elif spec.lower() == "auto":
        metrics = normalize_gene_column(
            read_table(gene_metrics_file), "gene metrics table"
        )
        required = {"fdr_phenotypes", "auroc_median", "nauprc_median"}
        missing = sorted(required - set(metrics.columns))
        if missing:
            raise ValueError(
                "Automatic panel selection requires metrics columns: "
                + ", ".join(missing)
            )
        if "status" in metrics.columns:
            metrics = metrics[metrics["status"].astype(str).eq("ok")].copy()
        for column in ("fdr_phenotypes", "auroc_median", "nauprc_median"):
            metrics[column] = pd.to_numeric(metrics[column], errors="coerce")
        fdr_rule = metrics["fdr_phenotypes"].gt(0)
        performance_rule = (
            metrics["auroc_median"].gt(0.8)
            & metrics["nauprc_median"].gt(10)
        )
        selected = metrics.loc[fdr_rule | performance_rule, ["Gene"]].copy()
        selected["Gene"] = selected["Gene"].astype(str).str.strip()
        requested = sorted(set(selected.loc[selected["Gene"].ne(""), "Gene"]))
    else:
        requested = read_gene_file(spec)

    selected_genes = [gene for gene in requested if gene in available_set]
    if not selected_genes:
        raise ValueError("No candidate-panel genes are present in the prediction table.")
    return selected_genes


def coerce_carrier(series):
    values = series.copy()
    mapped = values.astype(str).str.strip().str.lower().map(
        {"true": 1, "false": 0, "yes": 1, "no": 0}
    )
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric.where(mapped.isna(), mapped)
    invalid = numeric.notna() & ~numeric.isin([0, 1])
    return numeric.mask(invalid), invalid


def load_prediction_table(pred_file, gene_metrics_file, candidate_spec, top_fractions):
    predictions = read_table(pred_file)
    score_genes = {
        str(column)[:-6]
        for column in predictions.columns
        if str(column).endswith("_score")
        and f"{str(column)[:-6]}_carrier" in predictions.columns
    }
    if not score_genes:
        raise ValueError(
            "Prediction table must contain paired {GENE}_score and {GENE}_carrier columns."
        )
    panel = select_candidate_panel(
        candidate_spec, sorted(score_genes), gene_metrics_file
    )
    rows = []
    for gene in panel:
        scores = pd.to_numeric(predictions[f"{gene}_score"], errors="coerce")
        carriers, _ = coerce_carrier(predictions[f"{gene}_carrier"])
        valid = scores.notna() & carriers.notna()
        frame = pd.DataFrame({"score": scores[valid], "carrier": carriers[valid]})
        frame = frame.sort_values("score", ascending=False, kind="mergesort")
        n_total = int(len(frame))
        n_total_carriers = int(frame["carrier"].sum())
        for fraction in top_fractions:
            n_top = int(n_total * fraction)
            top = frame.iloc[:n_top]
            remaining = frame.iloc[n_top:]
            rows.append({
                "Gene": gene,
                "top_fraction": float(fraction),
                "top_percent": fraction_label(fraction),
                "percentile": percentile_label(fraction),
                "n_total": n_total,
                "n_total_carriers": n_total_carriers,
                "n_top": n_top,
                "n_top_carriers": int(top["carrier"].sum()),
                "n_top_noncarriers": int(n_top - top["carrier"].sum()),
                "n_remaining": int(len(remaining)),
                "n_remaining_carriers": int(remaining["carrier"].sum()),
                "n_remaining_noncarriers": int(
                    len(remaining) - remaining["carrier"].sum()
                ),
            })
    table = pd.DataFrame(rows)
    statistics = pd.DataFrame(
        [or_statistics(row) for _, row in table.iterrows()], index=table.index
    )
    table = pd.concat([table, statistics], axis=1)
    panel_order = {gene: index for index, gene in enumerate(panel)}
    table["_panel_order"] = table["Gene"].map(panel_order)
    table = table.sort_values(["_panel_order", "top_fraction"], kind="stable")
    return table.reset_index(drop=True), panel


def prepare_gene_or_details(table, panel):
    """Return the compact per-gene, per-threshold OR table."""
    panel_order = {gene: index for index, gene in enumerate(panel)}
    details = table.copy()
    details["_panel_order"] = details["Gene"].map(panel_order)
    details = details.sort_values(["_panel_order", "top_fraction"], kind="stable")
    preferred = [
        "Gene", "top_percent", "percentile", "n_total", "n_total_carriers",
        "n_top", "n_top_carriers", "or", "or_ci_low", "or_ci_upp", "p",
    ]
    return details[[column for column in preferred if column in details]].reset_index(
        drop=True
    )


def write_gene_file(path, genes):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(genes))
        handle.write("\n")


def extended_empirical_quantiles(values, probabilities=(0.25, 0.50, 0.75)):
    """Compute empirical quantiles while preserving zero and positive infinity."""
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    numeric = np.sort(numeric[~np.isnan(numeric)])
    if numeric.size == 0:
        return [np.nan] * len(probabilities)

    result = []
    for probability in probabilities:
        position = float(probability) * (numeric.size - 1)
        lower_index = int(math.floor(position))
        upper_index = int(math.ceil(position))
        lower_value = float(numeric[lower_index])
        upper_value = float(numeric[upper_index])
        if lower_index == upper_index or lower_value == upper_value:
            result.append(lower_value)
        elif np.isposinf(upper_value):
            result.append(np.inf)
        else:
            weight = position - lower_index
            result.append(lower_value + weight * (upper_value - lower_value))
    return result


def prepare_trend(table, top_fractions, finite_positive_only=False):
    available = pd.to_numeric(table["top_fraction"], errors="coerce").dropna()
    missing = [
        fraction for fraction in top_fractions
        if not np.isclose(available.to_numpy(dtype=float), fraction).any()
    ]
    if missing:
        labels = ", ".join(fraction_label(value) for value in missing)
        raise ValueError(
            "Requested top fraction(s) could not be computed from the prediction table: "
            f"{labels}."
        )
    rows = []
    for fraction in top_fractions:
        subset = table[np.isclose(table["top_fraction"], fraction)]
        odds_ratio = pd.to_numeric(subset["or"], errors="coerce")
        defined = odds_ratio.notna()
        finite_positive_mask = (
            defined & np.isfinite(odds_ratio) & odds_ratio.gt(0)
        )
        if finite_positive_only:
            included = finite_positive_mask
        else:
            # Keep defined boundary estimates: OR=0 and OR=infinity are
            # meaningful outcomes when a top or remaining group has no carriers.
            included = defined
        values = odds_ratio.loc[included].to_numpy(dtype=float)
        q1, median, q3 = extended_empirical_quantiles(values)
        rows.append(
            {
                "top_percent": fraction_label(fraction),
                "percentile": percentile_label(fraction),
                "n_genes": int(len(subset)),
                "or_median": float(median) if len(values) else np.nan,
                "or_q1": float(q1) if len(values) else np.nan,
                "or_q3": float(q3) if len(values) else np.nan,
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "top_percent", "percentile", "n_genes",
            "or_median", "or_q1", "or_q3",
        ],
    )


def lighten_color(color, white_fraction=0.70):
    rgb = np.asarray(mpl.colors.to_rgb(color), dtype=float)
    return tuple(rgb * (1.0 - white_fraction) + white_fraction)


def format_log_tick(value, _position):
    if value <= 0:
        return ""
    exponent = int(round(np.log10(value)))
    if np.isclose(value, 10.0**exponent):
        return f"{10.0**exponent:g}"
    return ""


def save_figure(figure, directory, stem):
    paths = []
    for extension in OUTPUT_FORMATS:
        path = os.path.join(directory, f"{stem}.{extension}")
        if extension == "png":
            figure.savefig(path, dpi=PLOT_DPI, bbox_inches="tight")
        else:
            figure.savefig(path, bbox_inches="tight")
        paths.append(path)
    plt.close(figure)
    return paths


def plot_median_trend(trend, figures_dir):
    data = trend.dropna(subset=["or_median", "or_q1", "or_q3"]).copy()
    if data.empty:
        raise ValueError("No OR values are available for the median OR trend figure.")

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.8,
        }
    )
    figure, axis = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    x = np.arange(len(data), dtype=float)
    median = data["or_median"].to_numpy(dtype=float)
    q1 = data["or_q1"].to_numpy(dtype=float)
    q3 = data["or_q3"].to_numpy(dtype=float)

    axis.plot(
        x,
        median,
        color=TREND_LINE_COLOR,
        linewidth=1.6,
        solid_capstyle="round",
        zorder=1,
    )
    for index, (value, lower, upper) in enumerate(zip(median, q1, q3)):
        axis.errorbar(
            index,
            value,
            yerr=[[value - lower], [upper - value]],
            fmt="D",
            markersize=5.0,
            markerfacecolor=lighten_color(TREND_OR_COLOR),
            markeredgecolor=TREND_OR_COLOR,
            markeredgewidth=1.1,
            ecolor=TREND_OR_COLOR,
            elinewidth=1.0,
            capsize=2.2,
            capthick=1.0,
            zorder=3,
        )

    axis.set_yscale("log", base=10)
    if (
        not np.isfinite(median).all()
        or not np.isfinite(q1).all()
        or not np.isfinite(q3).all()
        or (median <= 0).any()
        or (q1 <= 0).any()
        or (q3 <= 0).any()
    ):
        raise ValueError(
            "The requested raw OR summary contains a zero or infinity at a plotted "
            "quantile, which cannot be displayed on the log10 axis. The summary TSV "
            "was written successfully; rerun with --finite-positive-only for the "
            "finite-positive sensitivity plot."
        )
    lower_limit = min(1.0, max(0.1, float(np.nanmin(q1)) / 1.6))
    upper_limit = max(1.5, float(np.nanmax(q3)) * 1.6)
    axis.set_ylim(lower_limit, upper_limit)
    axis.set_xlim(-0.45, len(data) - 0.55)
    axis.set_xticks(x)
    axis.set_xticklabels(
        [percentile_tick_label(value) for value in data["percentile"]]
    )
    axis.set_xlabel("GPS percentile threshold")
    axis.set_ylabel("Median odds ratio (IQR)")

    axis.yaxis.set_major_locator(mticker.LogLocator(base=10, subs=(1.0,), numticks=6))
    axis.yaxis.set_major_formatter(mticker.FuncFormatter(format_log_tick))
    axis.yaxis.set_minor_locator(mticker.NullLocator())
    axis.yaxis.set_minor_formatter(mticker.NullFormatter())
    axis.tick_params(axis="y", which="minor", left=False, right=False)
    axis.axhline(
        1.0,
        color="#8A8A8A",
        linewidth=0.7,
        linestyle=(0, (3.0, 2.5)),
        zorder=0,
    )
    axis.grid(axis="y", which="major", color=GRID_COLOR, linewidth=0.5, zorder=0)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", which="major", direction="out", length=3.0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout(pad=0.35)
    return save_figure(
        figure,
        figures_dir,
        "fig_median_or_trend",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot the gene-level median OR trend for a fixed panel."
    )
    parser.add_argument("--outdir", default=None, help="Project output root.")
    parser.add_argument(
        "--pred-file", required=True,
        help="02 output containing spid, {GENE}_score and {GENE}_carrier columns.",
    )
    parser.add_argument(
        "--gene-metrics-file",
        required=True,
        help="Step 03 gene-level metrics table.",
    )
    parser.add_argument(
        "--candidate-genes",
        default="auto",
        help="Candidate panel: auto, all, or a one-gene-per-line file.",
    )
    parser.add_argument(
        "--top-percents",
        default=DEFAULT_TOP_PERCENTS,
        help=(
            "Comma-separated top fractions for the median OR trend; "
            "default: 0.2,0.1,0.05,0.01."
        ),
    )
    parser.add_argument(
        "--finite-positive-only",
        action="store_true",
        help=(
            "Restrict the median/IQR trend to finite positive OR values only. "
            "Without this switch, defined raw OR=0 and OR=infinity are retained."
        ),
    )
    args = parser.parse_args()

    project_output = os.path.abspath(args.outdir or OUTPUT_DIR)
    step_dir = os.path.join(project_output, "05_gene_level_or")
    output_dir = os.path.join(step_dir, "output")
    figures_dir = os.path.join(step_dir, "figures")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    print("Loading GPS predictions...")
    top_fractions = parse_top_percents(args.top_percents)
    table, panel = load_prediction_table(
        args.pred_file,
        args.gene_metrics_file,
        args.candidate_genes,
        top_fractions,
    )
    trend = prepare_trend(
        table,
        top_fractions,
        finite_positive_only=args.finite_positive_only,
    )
    details = prepare_gene_or_details(table, panel)
    trend_path = os.path.join(output_dir, "gene_level_or_median_trend.tsv")
    trend.to_csv(trend_path, sep="\t", index=False)
    details_path = os.path.join(output_dir, "gene_level_or_details.tsv")
    details.to_csv(details_path, sep="\t", index=False)
    candidate_path = os.path.join(output_dir, "selected_genes.txt")
    write_gene_file(candidate_path, panel)

    print("Plotting median OR trend...")
    figure_files = plot_median_trend(trend, figures_dir)

    print(f"Candidate gene list: {candidate_path}")
    print(f"Gene-level OR details: {details_path}")
    print(f"Gene-level OR figure: {figure_files[0]}")


if __name__ == "__main__":
    main()
