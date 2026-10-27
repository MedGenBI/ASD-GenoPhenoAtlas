import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


CODE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_ROOT not in sys.path:
    sys.path.insert(0, CODE_ROOT)

from src.config import (
    ALL_TRAITS,
    GENE_COL,
    MISSINGNESS_THRESHOLD,
    OUTPUT_DIR,
    PHENO_SAMPLE_COL,
)
from src.dataloader import (
    GENE_ASSOCIATION_COUNTS_FILENAME,
    normalize_phenotype_columns,
)


METRICS = ("auroc", "nauprc")
METRICS_OUTPUT_COLUMNS = [
    GENE_COL, "model", "n_valid", "n_carrier", "n_noncarriers",
    "prevalence", "auroc", "nauprc", "p_phenotypes", "fdr_phenotypes",
]
SUMMARY_OUTPUT_COLUMNS = [
    "metric", "n_genes_total", "n_genes_evaluable",
    "n_gps_above_diagonal", "n_general_severity_above_diagonal",
    "proportion_gps_above_diagonal", "paired_median_difference",
    "paired_max_difference",
]
DETAIL_OUTPUT_COLUMNS = [
    GENE_COL, "n_valid", "n_carrier", "prevalence", "metric",
    "gps_value", "general_severity_value", "gps_minus_general_severity",
]
REUSED_BOOTSTRAP_FIELDS = (
    "auroc_median", "auroc_ci_lower", "auroc_ci_upper",
    "nauprc_median", "nauprc_ci_lower", "nauprc_ci_upper",
)
FIGURE_DPI = 600
SIGNIFICANT_PHENOTYPE_COUNT_CMAP = LinearSegmentedColormap.from_list(
    "significant_phenotype_count",
    ["#709AE1", "#86C7C2", "#9BD29B", "#E3D375", "#F7E690"],
    N=256,
)


def infer_separator(path):
    return "," if str(path).lower().endswith(".csv") else "\t"


def read_table(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    return pd.read_csv(path, sep=infer_separator(path))


def write_table(frame, path, index=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    frame.to_csv(path, sep="\t", index=index)


def find_id_column(frame):
    names = {str(column).strip().lower(): column for column in frame.columns}
    for name in (
        PHENO_SAMPLE_COL.lower(),
        "sample_id",
        "sample id",
        "sample",
        "individual_id",
        "id",
    ):
        if name in names:
            return names[name]
    first = frame.columns[0]
    if str(first).lower().startswith("unnamed"):
        return first
    return first


def load_sample_table(path, normalize_traits=False):
    frame = read_table(path)
    id_column = find_id_column(frame)
    frame = frame.rename(columns={id_column: PHENO_SAMPLE_COL})
    if normalize_traits:
        frame = normalize_phenotype_columns(frame)
    frame[PHENO_SAMPLE_COL] = frame[PHENO_SAMPLE_COL].astype(str).str.strip()
    if frame[PHENO_SAMPLE_COL].duplicated().any():
        raise ValueError(f"Duplicate sample IDs in {path}")
    return frame.set_index(PHENO_SAMPLE_COL)


def extract_genes(predictions):
    genes = []
    for column in predictions.columns:
        if not str(column).endswith("_score"):
            continue
        gene = str(column)[:-6]
        if f"{gene}_carrier" in predictions.columns and gene not in genes:
            genes.append(gene)
    return genes


def coerce_carrier(series):
    values = series.copy()
    text = values.astype(str).str.strip().str.lower().map({
        "true": 1, "false": 0, "yes": 1, "no": 0,
    })
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric.where(text.isna(), text)
    invalid = numeric.notna() & ~numeric.isin([0, 1])
    return numeric.mask(invalid), invalid


def choose_binary_traits(phenotypes):
    configured_traits = [trait for trait in ALL_TRAITS if trait in phenotypes.columns]
    if not configured_traits:
        raise ValueError(
            "No configured phenotype traits were found in the phenotype file."
        )

    binary = {}
    invalid = []
    for column in configured_traits:
        values = pd.to_numeric(phenotypes[column], errors="coerce")
        non_numeric = phenotypes[column].notna() & values.isna()
        non_binary = values.notna() & ~values.isin([0, 1])
        if non_numeric.any() or non_binary.any():
            examples = phenotypes.loc[non_numeric | non_binary, column].astype(str).unique()[:4]
            invalid.append(f"{column}: {', '.join(examples)}")
            continue
        binary[column] = values
    if invalid:
        preview = "; ".join(invalid[:8])
        raise ValueError(
            "The phenotype file must contain only 0, 1, or missing values. "
            f"Invalid columns: {preview}"
        )
    if not binary:
        raise ValueError("No binary phenotype columns were found.")
    return pd.DataFrame(binary, index=phenotypes.index)


def severity_scores(binary):
    observed = binary.notna().sum(axis=1).astype(int)
    missing = binary.isna().sum(axis=1).astype(int)
    score = binary.eq(1).sum(axis=1).astype(float)
    score = score.mask(observed.eq(0))
    result = pd.DataFrame({
        "general_severity_score": score,
        "observed_phenotype_count": observed,
        "missing_phenotype_count": missing,
        "phenotype_missing_rate": missing / float(binary.shape[1]),
    }, index=binary.index)
    result.index.name = PHENO_SAMPLE_COL
    return result


def point_metrics(labels, scores):
    data = pd.concat([
        pd.to_numeric(labels, errors="coerce").rename("label"),
        pd.to_numeric(scores, errors="coerce").rename("score"),
    ], axis=1).dropna()
    n = len(data)
    carriers = int(data["label"].sum())
    prevalence = carriers / n if n else np.nan
    result = {
        "n_valid": int(n),
        "n_carrier": carriers,
        "n_noncarriers": int(n - carriers),
        "prevalence": prevalence,
        "auroc": np.nan,
        "nauprc": np.nan,
    }
    if n and data["label"].nunique() == 2:
        result["auroc"] = float(roc_auc_score(data["label"], data["score"]))
        average_precision = float(
            average_precision_score(data["label"], data["score"])
        )
        if prevalence > 0:
            result["nauprc"] = average_precision / prevalence
    return result, data


def empty_bootstrap_fields(source, reason):
    result = {field: np.nan for field in REUSED_BOOTSTRAP_FIELDS}
    result.update({
        "bootstrap_source": source,
        "bootstrap_validation_reason": reason,
    })
    return result


def load_bootstrap_rows(metrics_file, genes):
    if not metrics_file or not os.path.isfile(metrics_file):
        raise FileNotFoundError(f"03 gene-level metrics file not found: {metrics_file}")
    try:
        metrics = pd.read_csv(metrics_file, sep="\t")
    except Exception as error:
        raise ValueError(
            f"Unable to read the 03 gene-level metrics file: {metrics_file}"
        ) from error
    required = {
        GENE_COL, "n_samples", "n_carriers", "n_noncarriers", "prevalence",
        *REUSED_BOOTSTRAP_FIELDS,
    }
    missing = sorted(required - set(metrics.columns))
    if missing:
        raise ValueError(
            "The 03 gene-level metrics file is missing required columns: "
            + ", ".join(missing)
        )
    metrics["_gene_key"] = metrics[GENE_COL].astype(str)
    if metrics["_gene_key"].duplicated().any():
        duplicated = metrics.loc[metrics["_gene_key"].duplicated(), GENE_COL].astype(str).tolist()
        raise ValueError(
            "The 03 gene-level metrics file contains duplicated gene rows: "
            + ", ".join(duplicated[:10])
        )
    available = set(metrics["_gene_key"])
    missing_genes = [gene for gene in genes if str(gene) not in available]
    if missing_genes:
        raise ValueError(
            "The 03 gene-level metrics file is missing prediction genes: "
            + ", ".join(map(str, missing_genes[:10]))
        )
    rows = metrics.set_index("_gene_key").loc[[str(gene) for gene in genes]].copy()
    for field in ("n_samples", "n_carriers", "n_noncarriers", "prevalence"):
        values = pd.to_numeric(rows[field], errors="coerce")
        if values.isna().any():
            bad_genes = rows.index[values.isna()].tolist()
            raise ValueError(
                f"03 field '{field}' is not numeric for genes: "
                + ", ".join(map(str, bad_genes[:10]))
            )
    for field in REUSED_BOOTSTRAP_FIELDS:
        values = pd.to_numeric(rows[field], errors="coerce")
        finite = np.isfinite(values.to_numpy(dtype=float))
        if values.isna().any() or not finite.all():
            bad_genes = rows.index[values.isna() | ~finite].tolist()
            raise ValueError(
                f"03 bootstrap field '{field}' is not numeric for genes: "
                + ", ".join(map(str, bad_genes[:10]))
            )
    return rows


def validated_bootstrap_fields(row, gene, n_valid, n_carrier, n_noncarrier, prevalence):
    if int(row["n_samples"]) != int(n_valid):
        raise ValueError(
            f"03/current n_samples mismatch for {gene}: "
            f"03={row['n_samples']}, current={n_valid}"
        )
    if int(row["n_carriers"]) != int(n_carrier):
        raise ValueError(
            f"03/current n_carriers mismatch for {gene}: "
            f"03={row['n_carriers']}, current={n_carrier}"
        )
    if int(row["n_noncarriers"]) != int(n_noncarrier):
        raise ValueError(
            f"03/current n_noncarriers mismatch for {gene}: "
            f"03={row['n_noncarriers']}, current={n_noncarrier}"
        )
    if not np.isclose(float(row["prevalence"]), float(prevalence), atol=1e-12):
        raise ValueError(
            f"03/current prevalence mismatch for {gene}: "
            f"03={row['prevalence']}, current={prevalence}"
        )
    result = {field: row[field] for field in REUSED_BOOTSTRAP_FIELDS}
    result.update({
        "bootstrap_source": "03_reused",
        "bootstrap_validation_reason": "strictly_validated_against_current_paired_set",
    })
    return result


def add_association_counts(metrics, outdir):
    path = os.path.join(
        outdir, "01_data_preparation", "processed_data",
        GENE_ASSOCIATION_COUNTS_FILENAME,
    )
    if not os.path.isfile(path):
        return metrics, None
    counts = pd.read_csv(path)
    if GENE_COL not in counts.columns:
        counts = counts.rename(columns={counts.columns[0]: GENE_COL})
    keep = [
        column for column in (GENE_COL, "p_phenotypes", "fdr_phenotypes")
        if column in counts.columns
    ]
    if len(keep) == 1:
        return metrics, path
    return metrics.merge(
        counts[keep].drop_duplicates(GENE_COL), on=GENE_COL, how="left"
    ), path


def build_comparison(metrics):
    gps = metrics.loc[metrics["model"] == "gps"].set_index(GENE_COL)
    severity = metrics.loc[metrics["model"] == "general_severity"].set_index(GENE_COL)
    common = [gene for gene in gps.index if gene in severity.index]
    detail_rows = []
    for gene in common:
        gps_row = gps.loc[gene]
        severity_row = severity.loc[gene]
        for metric in METRICS:
            gps_value = gps_row[metric]
            severity_value = severity_row[metric]
            difference = gps_value - severity_value if pd.notna(gps_value) and pd.notna(severity_value) else np.nan
            detail_rows.append({
                GENE_COL: gene,
                "n_valid": gps_row["n_valid"],
                "n_carrier": gps_row["n_carrier"],
                "prevalence": gps_row["prevalence"],
                "metric": metric,
                "gps_value": gps_value,
                "general_severity_value": severity_value,
                "gps_minus_general_severity": difference,
            })
    details = pd.DataFrame(detail_rows)
    summary_rows = []
    for metric in METRICS:
        subset = details.loc[details["metric"] == metric]
        values = subset["gps_minus_general_severity"].dropna().astype(float)
        n = len(values)
        n_gps = int((values > 0).sum())
        n_severity = int((values < 0).sum())
        summary_rows.append({
            "metric": metric,
            "n_genes_total": int(len(subset)),
            "n_genes_evaluable": int(n),
            "n_gps_above_diagonal": n_gps,
            "n_general_severity_above_diagonal": n_severity,
            "proportion_gps_above_diagonal": n_gps / n if n else np.nan,
            "paired_median_difference": float(values.median()) if n else np.nan,
            "paired_max_difference": float(values.max()) if n else np.nan,
        })
    details = details[DETAIL_OUTPUT_COLUMNS]
    summary = pd.DataFrame(summary_rows, columns=SUMMARY_OUTPUT_COLUMNS)
    return details, summary


def parse_formats(value):
    formats = [item.strip().lower().lstrip(".") for item in value.split(",") if item.strip()]
    allowed = {"png", "pdf", "svg"}
    invalid = [item for item in formats if item not in allowed]
    if invalid:
        raise ValueError("Unsupported output format(s): " + ", ".join(invalid))
    return formats or ["png"]


def resolve_font():
    arial = [
        path for path in font_manager.findSystemFonts()
        if os.path.basename(path).lower().startswith(("arial", "arialmt"))
    ]
    if arial:
        font_manager.fontManager.addfont(arial[0])
        family = font_manager.FontProperties(fname=arial[0]).get_name()
    else:
        family = "DejaVu Sans"
    preferred_font = family if arial else "DejaVu Sans"
    plt.rcParams.update({
        "figure.dpi": 160,
        "savefig.dpi": FIGURE_DPI,
        "font.family": "sans-serif",
        "font.sans-serif": [preferred_font, "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "figure.titlesize": 7.0,
        "axes.titlesize": 7.0,
        "axes.labelsize": 7.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": "#D6DCE5",
        "grid.alpha": 0.62,
        "grid.linewidth": 0.45,
        "axes.edgecolor": "#4A4A4A",
        "axes.linewidth": 0.60,
        "xtick.color": "#4A4A4A",
        "ytick.color": "#4A4A4A",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def write_figure(fig, output_dir, name, formats):
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    for fmt in formats:
        path = os.path.join(output_dir, f"{name}.{fmt}")
        fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def to_wide(metrics):
    metrics = metrics.copy()
    metrics[GENE_COL] = metrics[GENE_COL].astype(str)
    if metrics.duplicated([GENE_COL, "model"]).any():
        raise ValueError("The gene-level metrics table contains duplicate Gene/model rows.")

    pieces = []
    for model in ("general_severity", "gps"):
        subset = metrics.loc[metrics["model"] == model].set_index(GENE_COL)
        if subset.empty:
            raise ValueError(f"Missing model rows: {model}")
        columns = [column for column in ("auroc", "nauprc", "n_carrier", "n_valid") if column in subset.columns]
        pieces.append(
            subset[columns].rename(columns={column: f"{model}_{column}" for column in columns})
        )

    wide = pieces[0].join(pieces[1], how="inner").reset_index()
    metadata_columns = [
        column for column in ("p_phenotypes", "fdr_phenotypes")
        if column in metrics.columns
    ]
    if metadata_columns:
        metadata = metrics[[GENE_COL] + metadata_columns].drop_duplicates(GENE_COL)
        wide = metadata.merge(wide, on=GENE_COL, how="right")
    return wide


def bubble_sizes(values, reference=None, min_size=8.0, max_size=80.0):
    current = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    ref = current if reference is None else pd.to_numeric(reference, errors="coerce").to_numpy(dtype=float)
    ref = ref[np.isfinite(ref)]
    if ref.size == 0:
        return np.full(current.shape, min_size, dtype=float)
    low, high = float(np.min(ref)), float(np.max(ref))
    if high <= low:
        return np.full(current.shape, (min_size + max_size) / 2.0, dtype=float)
    scaled = (current - low) / (high - low)
    scaled = np.clip(np.nan_to_num(scaled, nan=0.0), 0.0, 1.0)
    return min_size + (max_size - min_size) * np.sqrt(scaled)


def plot_gene_scatter(frame, output_dir, formats):
    specs = (
        ("auroc", "General severity AUROC", "GPS AUROC"),
        ("nauprc", "General severity nAUPRC", "GPS nAUPRC"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05), gridspec_kw={"wspace": 0.36})
    fig.subplots_adjust(left=0.085, right=0.855, bottom=0.19, top=0.965)
    color_column = next(
        (column for column in ("p_phenotypes", "nominal_phenotypes") if column in frame.columns),
        None,
    )
    size_column = "gps_n_carrier" if "gps_n_carrier" in frame.columns else None
    color_reference = pd.to_numeric(frame[color_column], errors="coerce").dropna() if color_column else None
    size_reference = pd.to_numeric(frame[size_column], errors="coerce").dropna() if size_column else None
    scatter_for_colorbar = None
    color_norm = None

    for axis, (metric, xlabel, ylabel) in zip(axes, specs):
        x = pd.to_numeric(frame[f"general_severity_{metric}"], errors="coerce")
        y = pd.to_numeric(frame[f"gps_{metric}"], errors="coerce")
        data = pd.DataFrame({"x": x, "y": y})
        if color_column:
            data["color"] = pd.to_numeric(frame[color_column], errors="coerce")
        if size_column:
            data["size"] = pd.to_numeric(frame[size_column], errors="coerce")
        data = data.replace([np.inf, -np.inf], np.nan).dropna()
        if data.empty:
            axis.set_visible(False)
            continue

        x = data["x"].astype(float)
        y = data["y"].astype(float)
        if metric == "nauprc":
            low = 0.0
            high = max(1.2, float(np.nanmax([x.max(), y.max()])) * 1.08)
            axis.set_xscale("symlog", linthresh=1.0, linscale=0.85, base=10)
            axis.set_yscale("symlog", linthresh=1.0, linscale=0.85, base=10)
            axis.axvline(1.0, color="#B7C1CB", linestyle=":", linewidth=0.7, zorder=0)
            axis.axhline(1.0, color="#B7C1CB", linestyle=":", linewidth=0.7, zorder=0)
        else:
            low = min(0.45, float(np.nanmin([x.min(), y.min()])))
            high = max(1.0, float(np.nanmax([x.max(), y.max()])))

        axis.plot([low, high], [low, high], color="#444444", linestyle="--", linewidth=1.0, zorder=1)
        point_sizes = bubble_sizes(data["size"], size_reference) if size_column else np.full(len(data), 42.0)
        if color_column:
            color_min = float(color_reference.min())
            color_max = float(color_reference.max())
            color_norm = PowerNorm(
                gamma=0.55,
                vmin=color_min,
                vmax=color_max if color_max > color_min else color_min + 1.0,
            )
            scatter_for_colorbar = axis.scatter(
                x, y, c=data["color"].astype(float),
                cmap=SIGNIFICANT_PHENOTYPE_COUNT_CMAP,
                norm=color_norm, s=point_sizes, alpha=0.82,
                edgecolor="#FFFFFF", linewidth=0.35, zorder=3,
            )
        else:
            axis.scatter(
                x, y, color="#557A9F", s=point_sizes, alpha=0.82,
                edgecolor="#FFFFFF", linewidth=0.35, zorder=3,
            )
        axis.set_xlim(low, high)
        axis.set_ylim(low, high)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.grid(color="#D6DCE5", linewidth=0.45, alpha=0.62)
        axis.set_axisbelow(True)
        axis.tick_params(width=0.6, length=2.5)

    if scatter_for_colorbar is not None:
        cbar_ax = fig.add_axes([0.895, 0.22, 0.009, 0.58])
        cbar_min = float(color_reference.min())
        cbar_max = float(color_reference.max())
        if cbar_max <= cbar_min:
            cbar_max = cbar_min + 1.0
        gradient = np.linspace(cbar_min, cbar_max, 256).reshape(-1, 1)
        cbar_ax.imshow(
            gradient, aspect="auto", origin="lower", extent=(0.0, 1.0, cbar_min, cbar_max),
            cmap=SIGNIFICANT_PHENOTYPE_COUNT_CMAP, norm=color_norm,
        )
        cbar_ax.set_xlim(0.0, 1.0)
        cbar_ax.set_ylim(cbar_min, cbar_max)
        cbar_ax.set_xticks([])
        tick_values = np.linspace(cbar_min, cbar_max, 6)
        if cbar_max - cbar_min >= 1:
            tick_values = np.unique(np.round(tick_values).astype(int))
        cbar_ax.set_yticks(tick_values)
        cbar_ax.set_ylabel("Nominal-significant phenotypes", fontsize=6.5, labelpad=4)
        cbar_ax.tick_params(labelsize=6.5, width=0.5, length=2)
        for spine in cbar_ax.spines.values():
            spine.set_linewidth(0.5)
            spine.set_color("#4A4A4A")
    return write_figure(fig, output_dir, "fig_gene_level_gps_vs_severity", formats)


def main():
    parser = argparse.ArgumentParser(
        description="Compare GPS with general binary severity and plot the results."
    )
    parser.add_argument("--outdir", default=None, help="Project output root.")
    parser.add_argument("--phenotype-file", required=True, help="Pre-binarized phenotype matrix.")
    parser.add_argument("--pred-file", required=True, help="GPS prediction TSV from step 02.")
    parser.add_argument("--gene-level-metrics-file", required=True, help="Gene-level metrics TSV from step 03.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated figure formats.")
    args = parser.parse_args()

    project_output = os.path.abspath(args.outdir or OUTPUT_DIR)
    step_dir = os.path.join(project_output, "04_gps_vs_severity")
    output_dir = os.path.join(step_dir, "output")
    figure_dir = os.path.join(step_dir, "figures")
    os.makedirs(output_dir, exist_ok=True)
    formats = parse_formats(args.formats)

    pred_path = os.path.abspath(args.pred_file)
    pheno_path = os.path.abspath(args.phenotype_file)
    metrics_path = os.path.abspath(args.gene_level_metrics_file)

    print(f"Loading predictions: {pred_path}")
    predictions = load_sample_table(pred_path)
    print(f"Loading binary phenotypes: {pheno_path}")
    phenotypes = load_sample_table(pheno_path, normalize_traits=True)
    genes = extract_genes(predictions)
    if not genes:
        raise ValueError("No score/carrier gene pairs found in the prediction file.")

    shared = predictions.index.intersection(phenotypes.index)
    if shared.empty:
        raise ValueError("Prediction and phenotype tables have no shared sample IDs.")
    predictions = predictions.loc[shared].copy()
    phenotypes = phenotypes.loc[shared].copy()
    binary = choose_binary_traits(phenotypes)
    missing_rate = binary.isna().mean(axis=1)
    keep = missing_rate <= MISSINGNESS_THRESHOLD
    predictions = predictions.loc[keep].copy()
    binary = binary.loc[keep].copy()
    if predictions.empty:
        raise ValueError("No samples remain after the 20% missingness filter.")

    severity = severity_scores(binary)
    bootstrap_rows = load_bootstrap_rows(metrics_path, genes)
    metrics_rows = []
    iterator = tqdm(genes, desc="Gene-level severity control", unit="gene") if tqdm else genes

    for gene in iterator:
        carrier, _ = coerce_carrier(predictions[f"{gene}_carrier"])
        gps = pd.to_numeric(predictions[f"{gene}_score"], errors="coerce")
        severity_score = severity["general_severity_score"]
        pairwise = carrier.notna() & gps.notna() & severity_score.notna()

        paired_labels = carrier.loc[pairwise]
        paired_gps = gps.loc[pairwise]
        paired_severity = severity_score.loc[pairwise]
        gps_point, _ = point_metrics(paired_labels, paired_gps)
        severity_point, _ = point_metrics(paired_labels, paired_severity)
        gps_bootstrap = validated_bootstrap_fields(
            bootstrap_rows.loc[str(gene)],
            gene,
            gps_point["n_valid"],
            gps_point["n_carrier"],
            gps_point["n_noncarriers"],
            gps_point["prevalence"],
        )

        common = {
            GENE_COL: gene,
            "evaluation_sample_mode": "paired",
            "n_total": len(predictions),
            "n_valid": gps_point["n_valid"],
            "n_carrier": gps_point["n_carrier"],
            "n_noncarriers": gps_point["n_noncarriers"],
            "prevalence": gps_point["prevalence"],
        }
        gps_record = {
            **common,
            "model": "gps",
            "auroc": gps_bootstrap["auroc_median"],
            "nauprc": gps_bootstrap["nauprc_median"],
            **gps_bootstrap,
        }
        severity_record = {
            **common,
            "model": "general_severity",
            "auroc": severity_point["auroc"],
            "nauprc": severity_point["nauprc"],
            **empty_bootstrap_fields("not_applicable", "severity_uses_point_estimate"),
        }
        metrics_rows.extend([gps_record, severity_record])

    metrics = pd.DataFrame(metrics_rows)
    metrics, _ = add_association_counts(metrics, project_output)
    details, summary = build_comparison(metrics)

    for column in ("p_phenotypes", "fdr_phenotypes"):
        if column not in metrics.columns:
            metrics[column] = np.nan
    metrics = metrics[METRICS_OUTPUT_COLUMNS]

    metrics_out = os.path.join(output_dir, "severity_control_gene_level_metrics.tsv")
    details_path = os.path.join(output_dir, "severity_control_paired_gene_differences.tsv")
    summary_path = os.path.join(output_dir, "severity_control_model_comparison_summary.tsv")

    write_table(metrics, metrics_out)
    write_table(details, details_path)
    write_table(summary, summary_path)

    print("Plotting GPS versus general severity...")
    resolve_font()
    plot_gene_scatter(to_wide(metrics), figure_dir, formats)
    print(f"Severity-control results and figures: {step_dir}")


if __name__ == "__main__":
    main()
