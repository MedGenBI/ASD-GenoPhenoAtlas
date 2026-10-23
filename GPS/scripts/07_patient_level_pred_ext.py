import argparse
import multiprocessing as mp
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product


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
    P_COL,
    PHENO_COL,
    PHENO_SAMPLE_COL,
    PLP_PATH,
    PLP_SAMPLE_COL,
    WGS_PHENO_PATH,
    WGS_PLP_PATH,
)
from src.dataloader import (
    PREPROCESSOR_FILENAME,
    load_training_data,
    normalize_phenotype_columns,
)
from src.model import GPSModel
from src.ranking import evaluate_rankings, read_gene_list, read_table


PARAM_GRID = list(product(
    np.round(np.arange(0.0, 1.0001, 0.05), 2),
    np.round(np.arange(0.0, 1.0001, 0.05), 2),
))
TUNING_STATE = {}
PREDICTION_STATE = {}


def panel_hash(genes):
    canonical = "\n".join(sorted(set(str(gene).strip() for gene in genes)))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def shared_trait_hash(derivation, external):
    shared = sorted(
        set(str(column) for column in derivation.columns)
        & set(str(column) for column in external.columns)
    )
    canonical = "\n".join(shared)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def save_hyperparameters(
    path,
    pen,
    conf,
    panel,
    aggregate,
    processed_dir,
    derivation,
    external,
):
    best = aggregate.iloc[0]
    payload = {
        "format": "gps_panel_hyperparameters_v1",
        "parameter_scope": "panel_global",
        "pen": float(pen),
        "conf": float(conf),
        "objective": "metric_val",
        "objective_definition": "adjusted AUPRC used by the public workflow",
        "parameter_grid": {
            "start": 0.0,
            "stop": 1.0,
            "step": 0.05,
        },
        "source_panel": {
            "gene_count": len(panel),
            "gene_hash": panel_hash(panel),
            "genes": list(panel),
        },
        "tuning_summary": {
            "mean_metric_val": float(best["mean_score"]),
            "n_valid_genes": int(best["n_valid_genes"]),
        },
        "derivation_processed_data": os.path.abspath(processed_dir),
        "shared_trait_columns_hash": shared_trait_hash(derivation, external),
    }
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_hyperparameters(path):
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Hyperparameter file must contain a JSON object.")
    scope = payload.get("parameter_scope")
    if scope is not None and scope != "panel_global":
        raise ValueError(
            "Hyperparameter file must contain panel_global parameters."
        )
    missing = [name for name in ("pen", "conf") if name not in payload]
    if missing:
        raise ValueError(
            "Hyperparameter file is missing: " + ", ".join(missing)
        )
    try:
        pen = float(payload["pen"])
        conf = float(payload["conf"])
    except (TypeError, ValueError) as error:
        raise ValueError("pen and conf must be numeric.") from error
    if not np.isfinite(pen) or not np.isfinite(conf):
        raise ValueError("pen and conf must be finite.")
    if not 0.0 <= pen <= 1.0 or not 0.0 <= conf <= 1.0:
        raise ValueError("pen and conf must be between 0 and 1.")
    return pen, conf, payload


def normalise_sample_table(frame, sample_column=PHENO_SAMPLE_COL):
    result = frame.copy()
    if sample_column not in result.columns:
        result = result.rename(columns={result.columns[0]: sample_column})
    result[sample_column] = result[sample_column].astype(str).str.strip()
    if result[sample_column].duplicated().any():
        raise ValueError("Phenotype table contains duplicate sample IDs.")
    return result.set_index(sample_column, drop=False)


def normalise_plp(frame):
    aliases = {
        "gene": GENE_COL,
        "gene_symbol": GENE_COL,
        "gene name": GENE_COL,
        "sample id": PHENO_SAMPLE_COL,
        "sample_id": PHENO_SAMPLE_COL,
        "spid": PHENO_SAMPLE_COL,
        "individual": PHENO_SAMPLE_COL,
        PLP_SAMPLE_COL.lower(): PHENO_SAMPLE_COL,
    }
    rename = {}
    for column in frame.columns:
        key = str(column).strip().lower()
        if key in aliases and aliases[key] not in frame.columns:
            rename[column] = aliases[key]
    result = frame.rename(columns=rename).copy()
    if not {GENE_COL, PHENO_SAMPLE_COL}.issubset(result.columns):
        raise ValueError("PLP table must contain gene and sample ID columns.")
    result[GENE_COL] = result[GENE_COL].astype(str).str.strip()
    result[PHENO_SAMPLE_COL] = result[PHENO_SAMPLE_COL].astype(str).str.strip()
    return result


def build_reference_tables(preprocessor):
    prevalence = pd.DataFrame(
        list(preprocessor.stats_prevalence.items()),
        columns=["phenotype", "prevalence"],
    )
    means = pd.DataFrame(
        list(preprocessor.stats_mean.items()),
        columns=["phenotype", "mean"],
    )
    return prevalence, means


def associated_phenotypes(association, gene, derivation, external):
    rows = association[association[GENE_COL].astype(str) == str(gene)].copy()
    p_values = pd.to_numeric(rows[P_COL], errors="coerce")
    available = set(derivation.columns) & set(external.columns)
    traits = rows.loc[
        p_values < ASSOCIATION_P_THRESHOLD, PHENO_COL
    ].dropna().astype(str)
    return sorted(set(trait for trait in traits if trait in available))


def adjusted_auprc(raw_auprc, prevalence):
    if not np.isfinite(raw_auprc) or not np.isfinite(prevalence):
        return np.nan
    if prevalence >= 1:
        return np.nan
    return float((raw_auprc - prevalence) / max(1e-12, 1.0 - prevalence))


def resolve_workers(value, n_tasks):
    value = int(value)
    if value < 1:
        raise ValueError("--workers must be a positive integer.")
    return max(1, min(value, n_tasks))


def process_context():
    if "fork" in mp.get_all_start_methods():
        return mp.get_context("fork")
    return mp.get_context()


def progress_iterator(futures, description):
    iterator = as_completed(futures)
    return (
        tqdm(iterator, total=len(futures), desc=description, unit="gene")
        if tqdm is not None
        else iterator
    )


def init_tuning_worker(
    derivation,
    association,
    derivation_plp,
    prevalence,
    means,
    trait_map,
):
    TUNING_STATE.clear()
    TUNING_STATE.update({
        "derivation": derivation,
        "association": association,
        "derivation_plp": derivation_plp,
        "prevalence": prevalence,
        "means": means,
        "trait_map": trait_map,
    })


def tune_one_gene(gene):
    context = TUNING_STATE
    derivation = context["derivation"]
    traits = context["trait_map"].get(gene, [])
    carriers = set(
        context["derivation_plp"].loc[
            context["derivation_plp"][GENE_COL] == gene, PHENO_SAMPLE_COL
        ]
    )
    if not traits:
        return [
            {
                GENE_COL: gene,
                "pen": pen,
                "conf": conf,
                "status": "skipped",
                "mean_score": np.nan,
                "n_samples": 0,
                "n_carriers": 0,
                "error": (
                    "No phenotype below the configured association threshold "
                    f"(P<{ASSOCIATION_P_THRESHOLD:g}) is available in both cohorts."
                ),
            }
            for pen, conf in PARAM_GRID
        ]

    valid = derivation[traits].isna().mean(axis=1) <= MISSINGNESS_THRESHOLD
    fit_data = derivation.loc[valid]
    y_true = fit_data[PHENO_SAMPLE_COL].astype(str).isin(carriers).astype(int).to_numpy()
    n_samples = len(y_true)
    n_carriers = int(y_true.sum())
    if n_samples == 0 or n_carriers == 0 or n_carriers == n_samples:
        return [
            {
                GENE_COL: gene,
                "pen": pen,
                "conf": conf,
                "status": "skipped",
                "mean_score": np.nan,
                "n_samples": n_samples,
                "n_carriers": n_carriers,
                "error": "Both carrier classes are required.",
            }
            for pen, conf in PARAM_GRID
        ]

    rows = []
    for pen, conf in PARAM_GRID:
        record = {
                GENE_COL: gene,
            "pen": pen,
            "conf": conf,
            "status": "failed",
            "mean_score": np.nan,
            "n_samples": n_samples,
            "n_carriers": n_carriers,
            "error": "",
        }
        try:
            model = GPSModel(
                gene,
                traits,
                context["association"],
                BINARY_TRAITS,
                CONTINUOUS_TRAITS,
                CONTINUOUS_TRAITS_DIRECTION,
            )
            model.fit(
                fit_data,
                pen,
                conf,
                context["prevalence"],
                context["means"],
            )
            scores = np.asarray(
                model.predict(
                    fit_data,
                    pen,
                    conf,
                    context["prevalence"],
                    context["means"],
                ),
                dtype=float,
            )
            raw_ap = float(average_precision_score(y_true, scores))
            record.update(
                status="ok",
                raw_auprc=raw_ap,
                prevalence=float(y_true.mean()),
                adjusted_auprc=adjusted_auprc(raw_ap, float(y_true.mean())),
                mean_score=adjusted_auprc(raw_ap, float(y_true.mean())),
            )
        except Exception as error:
            record["error"] = f"{type(error).__name__}: {str(error)[:400]}"
        rows.append(record)
    return rows


def tune_panel_global(
    genes,
    derivation,
    association,
    derivation_plp,
    prevalence,
    means,
    trait_map,
    workers,
):
    init_tuning_worker(
        derivation,
        association,
        derivation_plp,
        prevalence,
        means,
        trait_map,
    )
    results = []
    if workers == 1:
        iterator = (
            tqdm(genes, desc="Global hyperparameter optimization", unit="gene")
            if tqdm is not None else genes
        )
        for gene in iterator:
            results.extend(tune_one_gene(gene))
    else:
        kwargs = {
            "max_workers": workers,
            "initializer": init_tuning_worker,
            "initargs": (
                derivation,
                association,
                derivation_plp,
                prevalence,
                means,
                trait_map,
            ),
        }
        context = process_context()
        if context is not None:
            kwargs["mp_context"] = context
        with ProcessPoolExecutor(**kwargs) as executor:
            futures = {executor.submit(tune_one_gene, gene): gene for gene in genes}
            for future in progress_iterator(
                futures, "Global hyperparameter optimization"
            ):
                gene = futures[future]
                try:
                    results.extend(future.result())
                except Exception as error:
                    results.extend(
                        {
            GENE_COL: gene,
                            "pen": pen,
                            "conf": conf,
                            "status": "worker_failed",
                            "mean_score": np.nan,
                            "n_samples": 0,
                            "n_carriers": 0,
                            "error": f"{type(error).__name__}: {str(error)[:400]}",
                        }
                        for pen, conf in PARAM_GRID
                    )

    tuning = pd.DataFrame(results)
    valid = tuning[tuning["status"] == "ok"].copy()
    if valid.empty:
        raise ValueError(
            "Global hyperparameter optimization produced no valid results."
        )
    aggregate = (
        valid.groupby(["pen", "conf"], as_index=False)
        .agg(
            mean_score=("mean_score", "mean"),
            n_valid_genes=(GENE_COL, "nunique"),
        )
        .sort_values(["mean_score", "pen", "conf"], ascending=[False, True, True])
    )
    best = aggregate.iloc[0]
    return (
        float(best["pen"]),
        float(best["conf"]),
        tuning,
        aggregate,
    )


def init_prediction_worker(
    derivation,
    external,
    external_plp,
    derivation_plp,
    association,
    prevalence,
    means,
    trait_map,
    global_params,
):
    PREDICTION_STATE.clear()
    PREDICTION_STATE.update({
        "derivation": derivation,
        "external": external,
        "external_plp": external_plp,
        "derivation_plp": derivation_plp,
        "association": association,
        "prevalence": prevalence,
        "means": means,
        "trait_map": trait_map,
        "global_params": global_params,
    })


def predict_one_gene(gene):
    context = PREDICTION_STATE
    derivation = context["derivation"]
    external = context["external"]
    traits = context["trait_map"].get(gene, [])
    result = {
        GENE_COL: gene,
        "status": "failed",
        "n_derivation_samples": 0,
        "n_external_samples": len(external),
        "n_external_valid_scores": 0,
        "n_external_carriers": 0,
        "error": "",
        "scores": None,
        "carriers": None,
    }
    carriers = set(
        context["external_plp"].loc[
            context["external_plp"][GENE_COL] == gene, PHENO_SAMPLE_COL
        ]
    )
    carrier_series = pd.Series(
        external[PHENO_SAMPLE_COL].astype(str).isin(carriers).astype(int).to_numpy(),
        index=external.index,
        dtype="Int64",
    )
    result["carriers"] = carrier_series
    result["n_external_carriers"] = int(carrier_series.sum())
    if not traits:
        result["error"] = (
            "No phenotype below the configured association threshold "
            f"(P<{ASSOCIATION_P_THRESHOLD:g}) is available in both cohorts."
        )
        return result

    deriv_valid = derivation[traits].isna().mean(axis=1) <= MISSINGNESS_THRESHOLD
    external_valid = external[traits].isna().mean(axis=1) <= MISSINGNESS_THRESHOLD
    fit_data = derivation.loc[deriv_valid]
    prediction_data = external.loc[external_valid]
    y_true = fit_data[PHENO_SAMPLE_COL].astype(str).isin(
        set(
            context["derivation_plp"].loc[
                context["derivation_plp"][GENE_COL] == gene, PHENO_SAMPLE_COL
            ]
        )
    ).astype(int).to_numpy()
    if len(y_true) == 0 or y_true.sum() == 0 or y_true.sum() == len(y_true):
        result["error"] = "Derivation data does not contain both carrier classes."
        return result

    try:
        pen, conf = context["global_params"]
        model = GPSModel(
            gene,
            traits,
            context["association"],
            BINARY_TRAITS,
            CONTINUOUS_TRAITS,
            CONTINUOUS_TRAITS_DIRECTION,
        )
        model.fit(
            fit_data,
            pen,
            conf,
            context["prevalence"],
            context["means"],
        )
        score_values = model.predict(
            prediction_data,
            pen,
            conf,
            context["prevalence"],
            context["means"],
        )
        score_series = pd.Series(np.nan, index=external.index, dtype=float)
        score_series.loc[prediction_data.index] = np.asarray(
            score_values, dtype=float
        )
        result.update(
            status="ok",
            n_derivation_samples=len(fit_data),
            n_external_valid_scores=int(score_series.notna().sum()),
            scores=score_series,
        )
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {str(error)[:400]}"
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Predict GPS scores in an external validation cohort."
    )
    parser.add_argument(
        "--association-file",
        required=True,
        help="Association-results table generated by the GPA workflow.",
    )
    parser.add_argument("--outdir", default=None, help="Project output root.")
    parser.add_argument(
        "--candidate-genes",
        required=True,
        help="One-gene-per-line candidate panel file.",
    )
    parser.add_argument(
        "--external-pheno-file",
        default=None,
        help="External phenotype table; defaults to GPS_WGS_PHENO_PATH.",
    )
    parser.add_argument(
        "--external-plp-file",
        default=None,
        help="External PLP table; defaults to GPS_WGS_PLP_PATH.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of worker processes for tuning and prediction.",
    )
    parameter_group = parser.add_mutually_exclusive_group()
    parameter_group.add_argument(
        "--save-hyperparameters",
        default=None,
        help="Save the tuned panel-global parameters to a JSON file.",
    )
    parameter_group.add_argument(
        "--hyperparameters-file",
        default=None,
        help="Use fixed panel-global parameters from a JSON file and skip tuning.",
    )
    args = parser.parse_args()

    association_file = os.path.abspath(args.association_file)
    project_output = os.path.abspath(args.outdir or OUTPUT_DIR)
    step_dir = os.path.join(project_output, "07_patient_level_pred_ext")
    output_dir = os.path.join(step_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    external_pheno_file = os.path.abspath(args.external_pheno_file or WGS_PHENO_PATH)
    external_plp_file = os.path.abspath(args.external_plp_file or WGS_PLP_PATH)
    if not external_pheno_file or not os.path.isfile(external_pheno_file):
        raise FileNotFoundError(
            "Set GPS_WGS_PHENO_PATH or pass --external-pheno-file."
        )
    if not external_plp_file or not os.path.isfile(external_plp_file):
        raise FileNotFoundError(
            "Set GPS_WGS_PLP_PATH or pass --external-plp-file."
        )

    panel = read_gene_list(args.candidate_genes)
    if not panel:
        raise ValueError("The candidate gene panel is empty.")
    workers = resolve_workers(args.workers, len(panel))

    processed_dir = os.path.join(
        project_output, "01_data_preparation", "processed_data"
    )
    preprocessor_path = os.path.join(processed_dir, PREPROCESSOR_FILENAME)
    print("Loading derivation data...")
    print(f"Using derivation processed data: {processed_dir}")
    print(f"Using preprocessor: {preprocessor_path}")
    derivation, _, association, derivation_plp, preprocessor = load_training_data(
        association_file,
        PLP_PATH,
        processed_dir,
    )
    reference_values = getattr(preprocessor, "reference_values", None)
    if not isinstance(reference_values, dict):
        raise ValueError(
            "The preprocessor loaded from "
            f"{preprocessor_path} does not contain derivation ECDF references. "
            "Re-run scripts/01_data_preparation.py with the same --outdir "
            "before running this script."
        )
    missing_reference_traits = [
        trait
        for trait in CONTINUOUS_TRAITS
        if trait not in reference_values or len(reference_values[trait]) == 0
    ]
    if missing_reference_traits:
        raise ValueError(
            "The preprocessor loaded from "
            f"{preprocessor_path} is incomplete; missing ECDF references for: "
            + ", ".join(missing_reference_traits)
            + ". Re-run scripts/01_data_preparation.py with the same --outdir."
        )
    derivation = normalise_sample_table(derivation)
    derivation_plp = normalise_plp(derivation_plp)
    association = association.copy()
    association[GENE_COL] = association[GENE_COL].astype(str)
    association[PHENO_COL] = association[PHENO_COL].astype(str)

    print("Preprocessing external phenotypes with the frozen derivation reference...")
    external = normalise_sample_table(
        read_table(external_pheno_file),
    )
    external = normalize_phenotype_columns(external)
    external = preprocessor.transform(
        external,
        is_training=False,
        missing_threshold=MISSINGNESS_THRESHOLD,
    )
    external_plp = normalise_plp(read_table(external_plp_file))
    trait_map = {
        gene: associated_phenotypes(
            association,
            gene,
            derivation,
            external,
        )
        for gene in panel
    }
    prevalence, means = build_reference_tables(preprocessor)
    current_panel_hash = panel_hash(panel)
    current_trait_hash = shared_trait_hash(derivation, external)
    if args.hyperparameters_file:
        hyperparameters_path = os.path.abspath(args.hyperparameters_file)
        if not os.path.isfile(hyperparameters_path):
            raise FileNotFoundError(
                f"Hyperparameter file not found: {hyperparameters_path}"
            )
        pen, conf, metadata = load_hyperparameters(hyperparameters_path)
        source_panel = metadata.get("source_panel", {})
        source_hash = source_panel.get("gene_hash")
        if source_hash and source_hash != current_panel_hash:
            print(
                "Using parameters tuned on a different candidate panel "
                "(intentional fixed-parameter transfer)."
            )
        source_trait_hash = metadata.get("shared_trait_columns_hash")
        if source_trait_hash and source_trait_hash != current_trait_hash:
            print(
                "Warning: the shared trait-column set differs from the "
                "parameter source run."
            )
        print(f"Using fixed GPS hyperparameters: {hyperparameters_path}")
    else:
        print(
            f"Globally optimizing GPS hyperparameters across {len(panel)} genes..."
        )
        pen, conf, _, aggregate = tune_panel_global(
            panel,
            derivation,
            association,
            derivation_plp,
            prevalence,
            means,
            trait_map,
            workers,
        )
        if args.save_hyperparameters:
            hyperparameters_path = os.path.abspath(args.save_hyperparameters)
            save_hyperparameters(
                hyperparameters_path,
                pen,
                conf,
                panel,
                aggregate,
                processed_dir,
                derivation,
                external,
            )
            print(f"Saved GPS hyperparameters: {hyperparameters_path}")
    print(f"Predicting {len(panel)} genes for {len(external)} external samples...")
    init_prediction_worker(
        derivation,
        external,
        external_plp,
        derivation_plp,
        association,
        prevalence,
        means,
        trait_map,
        (pen, conf),
    )
    results = []
    if workers == 1:
        iterator = (
            tqdm(panel, desc="External prediction", unit="gene")
            if tqdm is not None else panel
        )
        for gene in iterator:
            result = predict_one_gene(gene)
            results.append(result)
    else:
        kwargs = {
            "max_workers": workers,
            "initializer": init_prediction_worker,
            "initargs": (
                derivation,
                external,
                external_plp,
                derivation_plp,
                association,
                prevalence,
                means,
                trait_map,
                (pen, conf),
            ),
        }
        context = process_context()
        if context is not None:
            kwargs["mp_context"] = context
        with ProcessPoolExecutor(**kwargs) as executor:
            futures = {executor.submit(predict_one_gene, gene): gene for gene in panel}
            for future in progress_iterator(futures, "External prediction"):
                gene = futures[future]
                try:
                    result = future.result()
                except Exception as error:
                    result = {
                        GENE_COL: gene,
                        "status": "worker_failed",
                        "n_derivation_samples": 0,
                        "n_external_samples": len(external),
                        "n_external_valid_scores": 0,
                        "n_external_carriers": 0,
                        "error": f"{type(error).__name__}: {str(error)[:400]}",
                        "scores": None,
                        "carriers": None,
                    }
                results.append(result)

    by_gene = {result[GENE_COL]: result for result in results}
    prediction = pd.DataFrame(index=external.index)
    prediction.index.name = PHENO_SAMPLE_COL
    for gene in panel:
        result = by_gene[gene]
        score = result.get("scores")
        if score is None:
            score = pd.Series(np.nan, index=external.index, dtype=float)
        carrier = result.get("carriers")
        if carrier is None:
            carrier = pd.Series(0, index=external.index, dtype="Int64")
        prediction[f"{gene}_score"] = score.reindex(external.index)
        prediction[f"{gene}_carrier"] = carrier.reindex(external.index)

    # Use the same panel-restricted ranking representation as the derivation step.
    print(f"Processing {len(panel)} panel genes across {len(external)} samples...")
    rank_samples = evaluate_rankings(prediction, panel, show_progress=True)
    rank_columns = [
        "Sample ID",
        "is_panel_carrier_sample",
        "n_rankable_panel_genes",
        "n_missing_panel_scores",
        "top1_gene",
        "top1_score",
        "true_carrier_genes",
        "true_carrier_scores",
        "best_true_rank",
    ]
    rank_samples = rank_samples[rank_samples["best_true_rank"].notna()].copy()
    rank_samples_path = os.path.join(output_dir, "patient_level_predication.tsv")
    rank_samples[rank_columns].to_csv(rank_samples_path, sep="\t", index=False)
    panel_path = os.path.join(output_dir, "patient_level_candidate_genes.txt")
    with open(panel_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(panel) + "\n")
    print(f"Patient-level ranking: {rank_samples_path}")


if __name__ == "__main__":
    main()
