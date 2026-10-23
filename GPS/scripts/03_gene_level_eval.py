import argparse
import hashlib
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


# Keep numerical libraries single-threaded inside each evaluation worker.
for name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(name, "1")

import numpy as np
import pandas as pd

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


CODE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_ROOT not in sys.path:
    sys.path.insert(0, CODE_ROOT)

from src.config import GENE_COL, OUTPUT_DIR, PHENO_SAMPLE_COL
from src.dataloader import (
    GENE_ASSOCIATION_COUNTS_FILENAME,
    read_table,
)
from src.metrics import bootstrap_auc
WORKER_STATE = {}

MAIN_COLUMNS = [
    GENE_COL, "n_samples", "n_carriers", "n_noncarriers", "prevalence",
    "auroc_median", "auroc_ci_lower", "auroc_ci_upper", "nauprc_median",
    "nauprc_ci_lower", "nauprc_ci_upper", "p_phenotypes", "fdr_phenotypes",
]

def resolve_workers(requested, n_tasks):
    if n_tasks < 1:
        return 1
    if requested is None or requested < 1:
        requested = os.cpu_count() or 1
    return max(1, min(int(requested), n_tasks))


def process_pool_context():
    if "fork" in mp.get_all_start_methods():
        return mp.get_context("fork")
    return None


def progress_futures(futures, description):
    iterator = as_completed(futures)
    if tqdm is not None:
        iterator = tqdm(iterator, total=len(futures), desc=description, unit="gene")
    return iterator


def deterministic_seed(base_seed, gene):
    payload = f"{int(base_seed)}:{gene}".encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    return int.from_bytes(digest, byteorder="little") % (2**32 - 1)


def coerce_carrier(series):
    values = series.copy()
    as_bool = values.astype(str).str.strip().str.lower().map({
        "true": 1, "false": 0, "yes": 1, "no": 0,
    })
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric.where(as_bool.isna(), as_bool)
    invalid = numeric.notna() & ~numeric.isin([0, 1])
    return numeric.mask(invalid), invalid


def base_record(gene):
    return {
        GENE_COL: gene,
        "n_samples": 0,
        "n_carriers": 0,
        "n_noncarriers": 0,
        "prevalence": np.nan,
        "auroc_median": np.nan,
        "auroc_ci_lower": np.nan,
        "auroc_ci_upper": np.nan,
        "nauprc_median": np.nan,
        "nauprc_ci_lower": np.nan,
        "nauprc_ci_upper": np.nan,
    }


def init_worker(predictions):
    global WORKER_STATE
    WORKER_STATE = {"predictions": predictions}


def evaluate_gene(gene, seed, n_bootstrap):
    predictions = WORKER_STATE["predictions"]
    score_col = f"{gene}_score"
    carrier_col = f"{gene}_carrier"
    record = base_record(gene)

    if score_col not in predictions.columns:
        return record
    if carrier_col not in predictions.columns:
        return record

    scores = pd.to_numeric(predictions[score_col], errors="coerce")
    carriers, _ = coerce_carrier(predictions[carrier_col])
    valid = scores.notna() & carriers.notna()
    frame = pd.DataFrame({
        "score": scores[valid].astype(float),
        "carrier": carriers[valid].astype(int),
    })
    record["n_samples"] = int(len(frame))
    if frame.empty:
        return record

    y_true = frame["carrier"].to_numpy(dtype=int)
    y_score = frame["score"].to_numpy(dtype=float)
    n_carriers = int(y_true.sum())
    record["n_carriers"] = n_carriers
    record["n_noncarriers"] = int(len(y_true) - n_carriers)
    record["prevalence"] = float(y_true.mean())
    if n_carriers == 0 or n_carriers == len(y_true):
        return record

    bootstrap = bootstrap_auc(
        y_true,
        y_score,
        n_bootstrap=n_bootstrap,
        random_state=deterministic_seed(seed, gene),
    )
    for field in (
        "auroc_median", "auroc_ci_lower", "auroc_ci_upper",
        "nauprc_median", "nauprc_ci_lower", "nauprc_ci_upper",
    ):
        record[field] = bootstrap[field]
    return record


def order_genes(genes, counts):
    order = {gene: index for index, gene in enumerate(counts)}
    return sorted(genes, key=lambda gene: (order.get(gene, len(order)), gene))


def build_summary(metrics):
    # A gene is evaluable only when the requested bootstrap median is available.
    good = metrics.loc[metrics["auroc_median"].notna()].copy()
    groups = {
        "all_evaluable_genes": np.ones(len(good), dtype=bool),
        "fdr_associated_genes": good["fdr_phenotypes"].fillna(0).gt(0).to_numpy(),
    }
    rows = []
    for group, mask in groups.items():
        subset = good.loc[mask]
        row = {
            "analysis_group": group,
            "n_genes": int(len(subset)),
            "auroc_median": subset["auroc_median"].median(),
            "nauprc_median": subset["nauprc_median"].median(),
        }
        rows.append(row)
    return pd.DataFrame(
        rows,
        columns=["analysis_group", "n_genes", "auroc_median", "nauprc_median"],
    )


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate GPS scores at the gene level."
    )
    parser.add_argument(
        "--pred-file", "-p", required=True,
        help="02 output containing spid, {GENE}_score and {GENE}_carrier columns.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of worker processes; the primary command uses 16.",
    )
    parser.add_argument(
        "--n-bootstrap", type=int, default=1000,
        help="Number of bootstrap resamples per gene; default: 1000.",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Base seed for deterministic per-gene bootstrap streams.",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Project output root; defaults to GPS_OUTPUT_DIR or ./output.",
    )
    args = parser.parse_args()
    if args.n_bootstrap < 1:
        raise ValueError("--n-bootstrap must be a positive integer.")
    project_output_dir = os.path.abspath(args.outdir or OUTPUT_DIR)
    step_dir = os.path.join(project_output_dir, "03_gene_level_eval")
    output_dir = os.path.join(step_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    counts_path = os.path.join(
        project_output_dir, "01_data_preparation", "processed_data",
        GENE_ASSOCIATION_COUNTS_FILENAME,
    )

    print(f"Project output: {project_output_dir}")
    print(f"Loading predictions: {args.pred_file}")
    predictions = pd.read_csv(args.pred_file, sep="\t")
    if PHENO_SAMPLE_COL not in predictions.columns:
        raise ValueError(
            f"Prediction file must contain a {PHENO_SAMPLE_COL} column."
        )
    if predictions[PHENO_SAMPLE_COL].duplicated().any():
        raise ValueError(
            f"Prediction file contains duplicated {PHENO_SAMPLE_COL} values."
        )
    predictions[PHENO_SAMPLE_COL] = predictions[PHENO_SAMPLE_COL].astype(str)
    predictions = predictions.set_index(PHENO_SAMPLE_COL)

    counts = read_table(counts_path)
    if GENE_COL not in counts.columns:
        counts = counts.rename(columns={counts.columns[0]: GENE_COL})
    counts[GENE_COL] = counts[GENE_COL].astype(str)
    counts_order = counts[GENE_COL].tolist()
    genes = [
        column[:-6] for column in predictions.columns
        if column.endswith("_score")
        and f"{column[:-6]}_carrier" in predictions.columns
    ]
    genes = order_genes(sorted(set(genes)), counts_order)
    if not genes:
        raise ValueError("No score/carrier gene pairs were found in the prediction file.")
    workers = resolve_workers(args.workers, len(genes))
    init_worker(predictions)

    print("Running gene-level evaluation...")
    results = []
    if workers == 1:
        iterator = tqdm(genes, desc="Gene-level evaluation", unit="gene") if tqdm else genes
        for gene in iterator:
            record = evaluate_gene(gene, args.seed, args.n_bootstrap)
            results.append(record)
    else:
        kwargs = {
            "max_workers": workers,
            "initializer": init_worker,
            "initargs": (predictions,),
        }
        context = process_pool_context()
        if context is not None:
            kwargs["mp_context"] = context
        with ProcessPoolExecutor(**kwargs) as executor:
            futures = {
                    executor.submit(
                        evaluate_gene,
                        gene,
                        args.seed,
                        args.n_bootstrap,
                    ): gene
                for gene in genes
            }
            for future in progress_futures(futures, "Gene-level evaluation"):
                gene = futures[future]
                try:
                    record = future.result()
                except Exception:
                    record = base_record(gene)
                results.append(record)

    metrics_full = pd.DataFrame(results)
    metrics_full = metrics_full.merge(
        counts[[GENE_COL, "p_phenotypes", "fdr_phenotypes"]],
        on=GENE_COL, how="left",
    )
    metrics_full = metrics_full.set_index(GENE_COL).loc[genes].reset_index()
    for column in MAIN_COLUMNS:
        if column not in metrics_full.columns:
            metrics_full[column] = np.nan
    metrics = metrics_full[MAIN_COLUMNS]

    metrics_path = os.path.join(output_dir, "gene_level_evaluation_metrics.tsv")
    summary_path = os.path.join(output_dir, "gene_level_evaluation_summary.tsv")
    metrics.to_csv(metrics_path, sep="\t", index=False)
    build_summary(metrics_full).to_csv(summary_path, sep="\t", index=False)

    print(f"Metrics: {metrics_path}")


if __name__ == "__main__":
    main()
