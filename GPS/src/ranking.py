import os

import numpy as np
import pandas as pd

from .config import PHENO_SAMPLE_COL

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


def infer_separator(path):
    return "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","


def read_table(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    return pd.read_csv(path, sep=infer_separator(path))


def read_prediction_table(path, sample_column=PHENO_SAMPLE_COL):
    frame = read_table(path)
    if sample_column not in frame.columns:
        frame = frame.rename(columns={frame.columns[0]: sample_column})
    frame[sample_column] = frame[sample_column].astype(str).str.strip()
    if frame[sample_column].duplicated().any():
        raise ValueError("Prediction table contains duplicate sample IDs.")
    return frame.set_index(sample_column)


def read_gene_list(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Gene panel file not found: {path}")

    genes = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            token = line.strip().split("\t")[0].split(",")[0].strip()
            if not token or token.startswith("#"):
                continue
            if token.lower() in {"gene", "gene_symbol", "genes"}:
                continue
            genes.append(token)
    return list(dict.fromkeys(genes))


def extract_genes(predictions):
    genes = []
    for column in predictions.columns:
        column = str(column)
        if column.endswith("_score"):
            gene = column[:-6]
            if f"{gene}_carrier" in predictions.columns and gene not in genes:
                genes.append(gene)
    return genes


def _numeric_matrix(frame, genes, suffix):
    columns = [f"{gene}{suffix}" for gene in genes]
    if not columns:
        return np.empty((len(frame), 0), dtype=float)
    return (
        frame.reindex(columns=columns)
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy(dtype=float)
    )


def _carrier_matrix(frame, genes):
    columns = [f"{gene}_carrier" for gene in genes]
    if not columns:
        return np.empty((len(frame), 0), dtype=float)

    raw = frame.reindex(columns=columns)
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    text = raw.astype("string").apply(
        lambda column: column.str.strip().str.lower().map(
            {"true": 1, "false": 0, "yes": 1, "no": 0}
        )
    )
    return numeric.where(text.isna(), text).to_numpy(dtype=float)


def _order_panel_scores(score_matrix, genes):
    """Return valid gene order by descending GPS score and ascending gene name."""
    n_samples, n_genes = score_matrix.shape
    if n_genes == 0:
        return (
            np.empty((n_samples, 0), dtype=int),
            np.empty((n_samples, 0), dtype=bool),
        )

    gene_order = np.argsort(np.asarray(genes, dtype=str), kind="stable")
    sorted_scores = score_matrix[:, gene_order]
    sort_key = np.where(np.isfinite(sorted_scores), -sorted_scores, np.inf)
    sorted_positions = np.argsort(sort_key, axis=1, kind="stable")
    ordered_indices = gene_order[sorted_positions]
    ordered_valid = np.take_along_axis(
        np.isfinite(score_matrix), ordered_indices, axis=1
    )
    return ordered_indices, ordered_valid


def evaluate_rankings(predictions, panel_genes, show_progress=False):
    """Create the sample-level ranking table consumed by scripts 06 and 07."""
    all_genes = extract_genes(predictions)
    available_genes = set(all_genes)
    panel = list(dict.fromkeys(gene for gene in panel_genes if gene in available_genes))
    if not panel:
        raise ValueError("No requested panel genes are available in the prediction table.")

    score_matrix = _numeric_matrix(predictions, panel, "_score")
    carrier_matrix = _carrier_matrix(predictions, panel)
    ordered_indices, ordered_valid = _order_panel_scores(score_matrix, panel)
    sample_ids = predictions.index.astype(str).to_numpy()

    row_iter = range(len(predictions))
    if show_progress and tqdm is not None:
        row_iter = tqdm(
            row_iter,
            total=len(predictions),
            desc="Ranking samples",
            unit="sample",
        )

    records = []
    for row_index in row_iter:
        carrier_mask = np.isfinite(carrier_matrix[row_index]) & (
            carrier_matrix[row_index] == 1
        )
        panel_true_indices = np.flatnonzero(carrier_mask)

        ranked_indices = ordered_indices[row_index][ordered_valid[row_index]]
        ranked_genes = [panel[index] for index in ranked_indices]
        ranked_scores = score_matrix[row_index, ranked_indices]
        position_by_gene = {
            gene: position
            for position, gene in enumerate(ranked_genes, start=1)
        }
        best_rank = min(
            (
                position_by_gene[panel[index]]
                for index in panel_true_indices
                if panel[index] in position_by_gene
            ),
            default=np.nan,
        )
        true_carrier_genes = [panel[index] for index in panel_true_indices]
        true_carrier_scores = [
            ""
            if not np.isfinite(score_matrix[row_index, index])
            else f"{float(score_matrix[row_index, index]):g}"
            for index in panel_true_indices
        ]

        records.append(
            {
                "Sample ID": sample_ids[row_index],
                "is_panel_carrier_sample": bool(panel_true_indices.size),
                "n_rankable_panel_genes": int(len(ranked_genes)),
                "n_missing_panel_scores": int(len(panel) - len(ranked_genes)),
                "top1_gene": ranked_genes[0] if ranked_genes else "",
                "top1_score": (
                    float(ranked_scores[0]) if len(ranked_scores) else np.nan
                ),
                "true_carrier_genes": ";".join(true_carrier_genes),
                "true_carrier_scores": ";".join(true_carrier_scores),
                "best_true_rank": (
                    float(best_rank) if np.isfinite(best_rank) else np.nan
                ),
            }
        )

    return pd.DataFrame(records)
