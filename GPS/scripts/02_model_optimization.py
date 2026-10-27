import argparse
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product


# Prevent nested BLAS threads when gene-level work is parallelized.
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
from sklearn.metrics import average_precision_score

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


CODE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_ROOT not in sys.path:
    sys.path.insert(0, CODE_ROOT)

from src.config import (
    ASSOCIATION_P_THRESHOLD,
    BINARY_TRAITS,
    CONTINUOUS_TRAITS,
    CONTINUOUS_TRAITS_DIRECTION,
    GENE_COL,
    MISSINGNESS_THRESHOLD,
    OUTPUT_DIR,
    PHENO_COL,
    PHENO_SAMPLE_COL,
    PLP_PATH,
    PLP_SAMPLE_COL,
    P_COL,
)
from src.dataloader import (
    load_training_data,
)
from src.model import GPSModel


METRIC_VALUE_NAME = "metric_val"

PENALTY_GRID = np.round(np.arange(0.0, 1.05, 0.05), 2)
CONFIDENCE_GRID = np.round(np.arange(0.0, 1.05, 0.05), 2)
PARAM_GRID = list(product(PENALTY_GRID, CONFIDENCE_GRID))

CACHE_COLUMNS = [
    "gene", "pen", "conf", "metric_val", "prevalence", "n_samples",
    "n_positive", "n_negative", "status", "error_type", "error_message",
]

PARAMETER_LOG_COLUMNS = [
    "gene", "best_pen", "best_conf", "cv_train_score",
    "cv_train_metric_val", "cv_train_metric_val_at_best", "tuning_gene_count",
    "prediction_n_samples", "prediction_n_carriers",
]

INNER_STATE = {}
PREDICTION_STATE = {}


def compute_metric_value(raw_auprc, prevalence):
    """Compute the fixed adjusted-AUPRC objective used by the public workflow."""
    if not np.isfinite(raw_auprc) or not np.isfinite(prevalence):
        return np.nan
    if prevalence < 0 or prevalence >= 1:
        return np.nan
    return float((raw_auprc - prevalence) / (1.0 - prevalence))


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
        iterator = tqdm(iterator, total=len(futures), desc=description, unit="task")
    return iterator


def build_gene_association_map(association, genes):
    work = association.copy()
    work["_gene"] = work[GENE_COL].astype(str)
    work["_p"] = pd.to_numeric(work[P_COL], errors="coerce")
    result = {}
    for gene in genes:
        traits = work.loc[
            (work["_gene"] == gene) & (work["_p"] < ASSOCIATION_P_THRESHOLD),
            PHENO_COL,
        ].dropna().astype(str).unique()
        result[gene] = sorted(traits.tolist())
    return result


def component_tables(preprocessor):
    prevalence = pd.DataFrame(
        list(preprocessor.stats_prevalence.items()),
        columns=["phenotype", "prevalence"],
    )
    means = pd.DataFrame(
        list(preprocessor.stats_mean.items()),
        columns=["phenotype", "mean"],
    )
    return prevalence, means


def empty_cache_record(gene, pen, conf, status, error_type="", error_message=""):
    return {
        "gene": gene,
        "pen": float(pen),
        "conf": float(conf),
        "metric_val": np.nan,
        "prevalence": np.nan,
        "n_samples": 0,
        "n_positive": 0,
        "n_negative": 0,
        "status": status,
        "error_type": error_type,
        "error_message": error_message,
    }


def init_inner_worker(df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map):
    global INNER_STATE
    INNER_STATE = {
        "df_pheno": df_pheno,
        "df_assoc": df_assoc,
        "df_plp": df_plp,
        "df_prev": df_prev,
        "df_mean": df_mean,
        "gene_assoc_map": gene_assoc_map,
    }


def compute_gene_inner_scores(gene):
    data = INNER_STATE
    phenotype = data["df_pheno"]
    association = data["df_assoc"]
    plp = data["df_plp"]
    prevalence = data["df_prev"]
    means = data["df_mean"]
    traits = data["gene_assoc_map"].get(gene, [])

    if not traits:
        return [
            empty_cache_record(
                gene, pen, conf, "no_associated_phenotypes",
                "NoAssociatedPhenotypes",
                "No nominally significant phenotype is available.",
            )
            for pen, conf in PARAM_GRID
        ]
    missing_traits = [trait for trait in traits if trait not in phenotype.columns]
    if missing_traits:
        message = "Missing phenotype columns: " + ",".join(missing_traits)
        return [
            empty_cache_record(gene, pen, conf, "missing_phenotype_columns",
                               "MissingColumn", message)
            for pen, conf in PARAM_GRID
        ]

    valid_mask = phenotype[traits].isna().mean(axis=1) <= MISSINGNESS_THRESHOLD
    valid_phenotype = phenotype.loc[valid_mask]
    carrier_ids = set(
        plp.loc[plp[GENE_COL].astype(str) == gene, PLP_SAMPLE_COL]
        .dropna().astype(str).unique()
    )
    if PHENO_SAMPLE_COL not in valid_phenotype.columns:
        message = f"Processed phenotype table must contain {PHENO_SAMPLE_COL}."
        return [
            empty_cache_record(gene, pen, conf, "invalid_input", "MissingColumn", message)
            for pen, conf in PARAM_GRID
        ]

    y_true = valid_phenotype[PHENO_SAMPLE_COL].astype(str).isin(carrier_ids).astype(int)
    n_samples = int(len(y_true))
    n_positive = int(y_true.sum())
    n_negative = n_samples - n_positive
    if n_samples == 0 or n_positive == 0 or n_negative == 0:
        status = "no_valid_rows" if n_samples == 0 else "single_class"
        error = "No valid rows." if n_samples == 0 else "Only one carrier class is present."
        return [
            {
                **empty_cache_record(gene, pen, conf, status, "SingleClass", error),
                "n_samples": n_samples,
                "n_positive": n_positive,
                "n_negative": n_negative,
            }
            for pen, conf in PARAM_GRID
        ]

    records = []
    for pen, conf in PARAM_GRID:
        try:
            model = GPSModel(
                gene,
                traits,
                association,
                BINARY_TRAITS,
                CONTINUOUS_TRAITS,
                CONTINUOUS_TRAITS_DIRECTION,
            )
            model.fit(valid_phenotype, pen, conf, prevalence, means)
            scores = model.predict(valid_phenotype, pen, conf, prevalence, means)
            scores = np.asarray(scores, dtype=float)
            raw_auprc = float(average_precision_score(y_true, scores))
            prevalence_value = float(y_true.mean())
            records.append({
                "gene": gene,
                "pen": float(pen),
                "conf": float(conf),
                "metric_val": compute_metric_value(raw_auprc, prevalence_value),
                "prevalence": prevalence_value,
                "n_samples": n_samples,
                "n_positive": n_positive,
                "n_negative": n_negative,
                "status": "ok",
                "error_type": "",
                "error_message": "",
            })
        except Exception as error:
            records.append({
                **empty_cache_record(
                    gene, pen, conf, "failed", type(error).__name__, str(error)[:500]
                ),
                "n_samples": n_samples,
                "n_positive": n_positive,
                "n_negative": n_negative,
            })
    return records


def compute_inner_cache(genes, workers, df_pheno, df_assoc, df_plp, df_prev, df_mean,
                        gene_assoc_map):
    init_inner_worker(df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map)
    results = []
    if workers == 1:
        iterator = tqdm(genes, desc="Inner gene scores", unit="gene") if tqdm else genes
        for gene in iterator:
            results.extend(compute_gene_inner_scores(gene))
    else:
        kwargs = {
            "max_workers": workers,
            "initializer": init_inner_worker,
            "initargs": (
                df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map
            ),
        }
        context = process_pool_context()
        if context is not None:
            kwargs["mp_context"] = context
        with ProcessPoolExecutor(**kwargs) as executor:
            futures = {executor.submit(compute_gene_inner_scores, gene): gene for gene in genes}
            for future in progress_futures(futures, "Inner gene scores"):
                gene = futures[future]
                try:
                    results.extend(future.result())
                except Exception as error:
                    results.extend(
                        empty_cache_record(
                            gene, pen, conf, "worker_failed",
                            type(error).__name__, str(error)[:500],
                        )
                        for pen, conf in PARAM_GRID
                    )
    return pd.DataFrame(results, columns=CACHE_COLUMNS).sort_values(
        ["gene", "pen", "conf"], kind="mergesort"
    ).reset_index(drop=True)


def select_best_parameters(cache, genes):
    """Apply leave-one-gene-out aggregation to the cached inner scores."""
    good = cache.loc[cache["status"] == "ok"].copy()
    selected = {}
    for test_gene in genes:
        train_genes = [gene for gene in genes if gene != test_gene]
        best_score = -np.inf
        best_pair = PARAM_GRID[0]
        best_rows = pd.DataFrame()
        for pen, conf in PARAM_GRID:
            rows = good.loc[
                good["gene"].isin(train_genes)
                & (good["pen"] == float(pen))
                & (good["conf"] == float(conf))
            ]
            score = rows[METRIC_VALUE_NAME].mean()
            score = float(score) if np.isfinite(score) else 0.0
            if score > best_score:
                best_score = score
                best_pair = (float(pen), float(conf))
                best_rows = rows

        train_total = len(train_genes)
        metric_value = (
            float(best_rows[METRIC_VALUE_NAME].mean())
            if not best_rows.empty else np.nan
        )
        selected[test_gene] = {
            "gene": test_gene,
            "best_pen": best_pair[0],
            "best_conf": best_pair[1],
            "cv_train_score": best_score,
            "cv_train_metric_val": metric_value,
            "cv_train_metric_val_at_best": metric_value,
            "tuning_gene_count": train_total,
            "prediction_n_samples": np.nan,
            "prediction_n_carriers": np.nan,
        }
    return selected


def init_prediction_worker(df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map):
    global PREDICTION_STATE
    PREDICTION_STATE = {
        "df_pheno": df_pheno,
        "df_assoc": df_assoc,
        "df_plp": df_plp,
        "df_prev": df_prev,
        "df_mean": df_mean,
        "gene_assoc_map": gene_assoc_map,
    }


def predict_one_gene(gene, best_info):
    data = PREDICTION_STATE
    phenotype = data["df_pheno"]
    association = data["df_assoc"]
    plp = data["df_plp"]
    prevalence = data["df_prev"]
    means = data["df_mean"]
    traits = data["gene_assoc_map"].get(gene, [])
    result = dict(best_info)
    result.update({
        "prediction_status": "pending",
        "prediction_error_type": "",
        "prediction_error_message": "",
        "pred_df": None,
    })
    try:
        if not traits:
            raise ValueError("No nominally significant associated phenotype.")
        valid_mask = phenotype[traits].isna().mean(axis=1) <= MISSINGNESS_THRESHOLD
        valid_phenotype = phenotype.loc[valid_mask]
        model = GPSModel(
            gene,
            traits,
            association,
            BINARY_TRAITS,
            CONTINUOUS_TRAITS,
            CONTINUOUS_TRAITS_DIRECTION,
        )
        model.fit(
            valid_phenotype,
            best_info["best_pen"],
            best_info["best_conf"],
            prevalence,
            means,
        )
        scores = model.predict(
            valid_phenotype,
            best_info["best_pen"],
            best_info["best_conf"],
            prevalence,
            means,
        )
        carriers = set(
            plp.loc[plp[GENE_COL].astype(str) == gene, PLP_SAMPLE_COL]
            .dropna().astype(str).unique()
        )
        result["pred_df"] = pd.DataFrame({
            PHENO_SAMPLE_COL: valid_phenotype[PHENO_SAMPLE_COL].astype(str),
            f"{gene}_score": np.asarray(scores, dtype=float),
            f"{gene}_carrier": valid_phenotype[PHENO_SAMPLE_COL].astype(str).isin(carriers).astype(int),
        })
        result["prediction_status"] = "ok"
        result["prediction_n_samples"] = int(len(valid_phenotype))
        result["prediction_n_carriers"] = int(result["pred_df"][f"{gene}_carrier"].sum())
    except Exception as error:
        result["prediction_status"] = "failed"
        result["prediction_error_type"] = type(error).__name__
        result["prediction_error_message"] = str(error)[:500]
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run the fixed primary GPS LOG-CV workflow."
    )
    parser.add_argument(
        "--association-file",
        required=True,
        help="Association-results table generated by the GPA workflow.",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of worker processes; the primary command uses 16.",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Project output root; defaults to GPS_OUTPUT_DIR or ./output.",
    )
    args = parser.parse_args()

    association_file = os.path.abspath(args.association_file)
    if not PLP_PATH:
        raise ValueError("Set GPS_PLP_PATH before running this script.")

    project_output_dir = os.path.abspath(args.outdir or OUTPUT_DIR)
    input_processed_data_dir = os.path.join(
        project_output_dir, "01_data_preparation", "processed_data"
    )
    step_dir = os.path.join(project_output_dir, "02_model_optimization")
    output_dir = os.path.join(step_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Project output: {project_output_dir}")
    print(f"Input processed data: {input_processed_data_dir}")
    print(f"02 output: {output_dir}")

    df_pheno, df_gene_counts, df_assoc, df_plp, preprocessor = load_training_data(
        assoc_path=association_file,
        plp_path=PLP_PATH,
        processed_data_dir=input_processed_data_dir,
    )
    df_prev, df_mean = component_tables(preprocessor)
    df_gene_counts.index = df_gene_counts.index.astype(str)
    genes = sorted(
        df_gene_counts.loc[df_gene_counts["p_phenotypes"] > 0].index.tolist()
    )
    if not genes:
        raise ValueError("No genes with nominally significant associations were found.")
    gene_assoc_map = build_gene_association_map(df_assoc, genes)
    workers = resolve_workers(args.workers, len(genes))

    print(f"Computing {len(genes)} x {len(PARAM_GRID)} inner scores...")
    cache = compute_inner_cache(
        genes, workers, df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map
    )
    cache_path = os.path.join(output_dir, "gene_parameter_inner_scores.tsv")
    cache.to_csv(cache_path, sep="\t", index=False)
    print("Running LOG-CV...")
    best_by_gene = select_best_parameters(cache, genes)

    init_prediction_worker(
        df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map
    )
    results = []
    if workers == 1:
        iterator = tqdm(genes, desc="Outer gene predictions", unit="gene") if tqdm else genes
        for gene in iterator:
            results.append(predict_one_gene(gene, best_by_gene[gene]))
    else:
        kwargs = {
            "max_workers": workers,
            "initializer": init_prediction_worker,
            "initargs": (
                df_pheno, df_assoc, df_plp, df_prev, df_mean, gene_assoc_map
            ),
        }
        context = process_pool_context()
        if context is not None:
            kwargs["mp_context"] = context
        with ProcessPoolExecutor(**kwargs) as executor:
            futures = {
                executor.submit(predict_one_gene, gene, best_by_gene[gene]): gene
                for gene in genes
            }
            for future in progress_futures(futures, "Outer gene predictions"):
                gene = futures[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    failed = dict(best_by_gene[gene])
                    failed.update({
                        "prediction_status": "failed",
                        "prediction_error_type": type(error).__name__,
                        "prediction_error_message": str(error)[:500],
                        "pred_df": None,
                    })
                    results.append(failed)

    results.sort(key=lambda item: item["gene"])
    parameter_log_full = pd.DataFrame([
        {key: value for key, value in result.items() if key != "pred_df"}
        for result in results
    ])
    parameter_log = parameter_log_full.reindex(columns=PARAMETER_LOG_COLUMNS)
    parameter_path = os.path.join(output_dir, "logcv_best_params.csv")
    parameter_log.to_csv(parameter_path, index=False)

    prediction_tables = [
        result["pred_df"].set_index(PHENO_SAMPLE_COL)
        for result in results if result.get("pred_df") is not None
    ]
    prediction_path = os.path.join(output_dir, "gps_predictions_logcv.tsv")
    if prediction_tables:
        pd.concat(prediction_tables, axis=1).sort_index().to_csv(
            prediction_path, sep="\t"
        )
        print(f"Predictions: {prediction_path}")
    else:
        print("No valid gene prediction was produced.")


if __name__ == "__main__":
    main()
