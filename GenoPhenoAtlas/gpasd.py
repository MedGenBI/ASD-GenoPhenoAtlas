from __future__ import annotations

import argparse
import os
import time
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import scipy.stats as stats

ID_COLUMN = "spid"
PHENOTYPE_SEPARATOR = ","
GENOTYPE_SEPARATOR = "\t"
COVARIATE_SEPARATOR = ","
MIN_CARRIERS = 5
MIN_CASES = 5
PERMUTATION_COUNT = 1000
PERMUTATION_SEED = 42
TARGET_FDR = 0.05
DEFAULT_TRAIT_METADATA = "Data/trait_metadata.csv"
TRAIT_METADATA_COLUMNS = (
    "trait_id",
    "trait_label",
    "trait_type",
    "clinical_domain",
    "effect_direction",
)
TRAIT_TYPES = {"binary", "quantitative"}
EFFECT_DIRECTIONS = {"higher", "lower"}

ASSOCIATION_COLUMNS = [
    "Trait",
    "Gene",
    "Beta",
    "SE",
    "OR",
    "95% CI Lower",
    "95% CI Upper",
    "P",
    "N Carriers",
    "Trait category",
    "Analysis method",
    "N Cases",
    "Affected cases",
    "Unaffected cases",
    "Note",
]


def progress(iterable, total: int, label: str):
    step = max(1, total // 20)
    for index, item in enumerate(iterable, 1):
        if index == 1 or index == total or index % step == 0:
            print(f"{label}: {index}/{total}", flush=True)
        yield item


def ensure_parent(path: str | Path) -> None:
    parent = Path(path).parent
    if str(parent) != ".":
        parent.mkdir(parents=True, exist_ok=True)


def read_matrix(path: str, sep: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=sep, index_col=ID_COLUMN)
    if frame.index.has_duplicates:
        raise ValueError(f"Duplicate sample IDs found in {path}.")
    return frame


def load_trait_metadata(path: str = DEFAULT_TRAIT_METADATA) -> pd.DataFrame:
    metadata_path = Path(path)
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Trait metadata file not found: {metadata_path}")

    metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
    missing_columns = set(TRAIT_METADATA_COLUMNS).difference(metadata.columns)
    if missing_columns:
        raise ValueError(
            "Trait metadata is missing required columns: "
            + ", ".join(sorted(missing_columns))
            + "."
        )

    metadata = metadata[list(TRAIT_METADATA_COLUMNS)].copy()
    for column in metadata.columns:
        metadata[column] = metadata[column].str.strip()

    required_values = ["trait_id", "trait_label", "trait_type", "clinical_domain"]
    missing_values = {
        column
        for column in required_values
        if metadata[column].eq("").any()
    }
    if missing_values:
        raise ValueError(
            "Trait metadata contains missing values in: "
            + ", ".join(sorted(missing_values))
            + "."
        )
    if metadata["trait_id"].duplicated().any():
        raise ValueError("Trait metadata contains duplicate trait_id values.")

    invalid_types = sorted(set(metadata["trait_type"]) - TRAIT_TYPES)
    if invalid_types:
        raise ValueError(
            "Trait metadata contains invalid trait_type values: "
            + ", ".join(invalid_types)
            + "."
        )

    quantitative = metadata["trait_type"].eq("quantitative")
    invalid_directions = sorted(
        set(metadata.loc[quantitative, "effect_direction"]) - EFFECT_DIRECTIONS
    )
    if invalid_directions:
        raise ValueError(
            "Quantitative traits must use effect_direction 'higher' or 'lower'."
        )

    return metadata


def trait_metadata_maps(
    metadata: pd.DataFrame,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    indexed = metadata.set_index("trait_id")
    return (
        indexed["trait_type"].to_dict(),
        indexed["trait_label"].to_dict(),
        indexed["clinical_domain"].to_dict(),
        indexed["effect_direction"].to_dict(),
    )


def validate_traits(traits: object, trait_types: dict[str, str]) -> None:
    unknown = sorted(set(traits) - set(trait_types))
    if unknown:
        raise ValueError(
            "Trait metadata is missing definitions for: " + ", ".join(unknown)
        )


def prepare_covariates(frame: pd.DataFrame) -> pd.DataFrame:
    categorical = frame.select_dtypes(include=["object", "category", "string"]).columns.tolist()
    if categorical:
        frame = pd.get_dummies(frame, columns=categorical, drop_first=True, dtype=float)

    frame = frame.astype(float)
    zero_variance = frame.columns[frame.var() == 0].tolist()
    if zero_variance:
        frame = frame.drop(columns=zero_variance)

    scale_columns = [column for column in frame.columns if frame[column].nunique() > 2]
    if scale_columns:
        frame[scale_columns] = (
            frame[scale_columns] - frame[scale_columns].mean()
        ) / frame[scale_columns].std()
    return frame


def prepare_observed_inputs(
    phenotype: pd.DataFrame,
    genotype: pd.DataFrame,
    covariate: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    common = phenotype.index.intersection(genotype.index).intersection(covariate.index)
    phenotype = phenotype.loc[common].copy()
    genotype = genotype.loc[common].copy()
    covariate = covariate.loc[common].copy()

    complete_covariates = ~covariate.isnull().any(axis=1)
    phenotype = phenotype.loc[complete_covariates].astype(float)
    genotype = genotype.loc[complete_covariates].astype(float)
    covariate = prepare_covariates(covariate.loc[complete_covariates])
    return phenotype, genotype, covariate


def prepare_permutation_inputs(
    phenotype: pd.DataFrame,
    covariate: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = phenotype.index.intersection(covariate.index)
    phenotype = phenotype.loc[common].copy()
    covariate = covariate.loc[common].copy()
    complete_covariates = ~covariate.isnull().any(axis=1)
    phenotype = phenotype.loc[complete_covariates].astype(float)
    covariate = prepare_covariates(covariate.loc[complete_covariates])
    return phenotype, covariate


def inverse_normal_transform(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    transformed = np.full(values.shape, np.nan)
    valid = ~np.isnan(values)
    ranks = stats.rankdata(values[valid], method="average")
    n = len(ranks)
    if n == 0:
        return transformed
    c = 3 / 8
    probabilities = (ranks - c) / (n - 2 * c + 1)
    transformed[valid] = stats.norm.ppf(probabilities)
    return transformed


def fit_ols(design: np.ndarray, outcome: np.ndarray) -> tuple[float, float, float]:
    n, k = design.shape
    if n <= k:
        return np.nan, np.nan, np.nan
    try:
        xtx_inverse = np.linalg.inv(design.T @ design)
        beta = xtx_inverse @ design.T @ outcome
        residual = outcome - design @ beta
        residual_variance = np.sum(residual**2) / (n - k)
        standard_errors = np.sqrt(np.diag(xtx_inverse) * residual_variance)
        effect = beta[1]
        standard_error = standard_errors[1]
        statistic = effect / standard_error
        p_value = 2 * (1 - stats.t.cdf(np.abs(statistic), df=n - k))
        return effect, standard_error, p_value
    except (FloatingPointError, np.linalg.LinAlgError, ValueError):
        return np.nan, np.nan, np.nan


def fit_firth(design: np.ndarray, outcome: np.ndarray) -> tuple[float, float, float, str]:
    try:
        from firthmodels import FirthLogisticRegression
    except ImportError as error:
        raise RuntimeError("The firthmodels package is required for binary traits.") from error

    note = "Success"
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        model = FirthLogisticRegression(fit_intercept=False)
        model.fit(design, outcome).lrt()
        for warning in captured:
            if "Maximum number of iterations reached" in str(warning.message):
                note = "Firth non-convergence warning"
                break

    coefficients = np.atleast_1d(np.squeeze(model.coef_))
    standard_errors = np.atleast_1d(np.squeeze(model.lrt_bse_))
    p_values = np.atleast_1d(np.squeeze(model.lrt_pvalues_))
    return coefficients[1], standard_errors[1], p_values[1], note


def fit_pair(
    trait: str,
    gene: str,
    genotype: np.ndarray,
    outcome: np.ndarray,
    covariates: np.ndarray,
    continuous: bool,
) -> dict[str, object] | None:
    valid_genotype = ~np.isnan(genotype)
    genotype = genotype[valid_genotype]
    outcome = outcome[valid_genotype]
    covariates = covariates[valid_genotype]

    n_carriers = int(np.sum(genotype == 1))
    if n_carriers < MIN_CARRIERS:
        return None

    n_samples = len(outcome)
    result: dict[str, object] = {
        "Trait": trait,
        "Gene": gene,
        "Beta": np.nan,
        "SE": np.nan,
        "OR": np.nan,
        "95% CI Lower": np.nan,
        "95% CI Upper": np.nan,
        "P": np.nan,
        "N Carriers": n_carriers,
        "Trait category": "quantitative" if continuous else "binary",
        "Analysis method": "linear" if continuous else "firth logistic",
        "N Cases": n_samples,
        "Affected cases": np.nan,
        "Unaffected cases": np.nan,
        "Note": "Success",
    }

    if n_samples < 50:
        result["Note"] = "Skipped: total sample size < 50"
        return result

    design = np.column_stack((np.ones(n_samples), genotype, covariates))

    if continuous:
        beta, standard_error, p_value = fit_ols(design, outcome)
        result["Beta"] = beta
        result["SE"] = standard_error
        result["P"] = p_value
        if np.isnan(p_value):
            result["Note"] = "OLS fitting failed"
        return result

    binary_outcome = outcome.astype(int)
    affected = int(np.sum(binary_outcome == 1))
    unaffected = int(np.sum(binary_outcome == 0))
    result["Affected cases"] = affected
    result["Unaffected cases"] = unaffected
    if affected < MIN_CASES or unaffected < MIN_CASES:
        result["Note"] = "Skipped: fewer than five affected or unaffected individuals"
        return result

    carrier_cases = int(np.sum((genotype > 0) & (binary_outcome == 1)))
    carrier_controls = int(np.sum((genotype > 0) & (binary_outcome == 0)))

    try:
        beta, standard_error, p_value, note = fit_firth(design, binary_outcome)
        result["Beta"] = beta
        result["SE"] = standard_error
        result["P"] = p_value
        result["Note"] = note
        result["OR"] = np.exp(beta)
        result["95% CI Lower"] = np.exp(beta - 1.96 * standard_error)
        result["95% CI Upper"] = np.exp(beta + 1.96 * standard_error)
    except Exception as error:
        result["Note"] = f"Firth fitting failed: {error}"

    if carrier_cases == 0 or carrier_controls == 0:
        result["Note"] = f"{result['Note']}; complete separation handled by Firth"
    return result


def phenotype_arrays(
    phenotype: np.ndarray,
    genotype: np.ndarray,
    covariates: np.ndarray,
    continuous: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid_phenotype = ~np.isnan(phenotype)
    outcome = phenotype[valid_phenotype]
    genotype = genotype[valid_phenotype]
    covariates = covariates[valid_phenotype]
    if continuous:
        outcome = inverse_normal_transform(outcome)
    else:
        outcome = outcome.astype(int)
    return outcome, genotype, covariates


def process_phenotype(task: tuple[object, ...]) -> list[dict[str, object]]:
    trait, phenotype, genotype, covariates, trait_types = task
    continuous = trait_types[trait] == "quantitative"
    outcome, genotype_values, covariate_values = phenotype_arrays(
        phenotype, genotype.to_numpy(), covariates, continuous
    )
    results = []
    for column_index, gene in enumerate(genotype.columns):
        result = fit_pair(
            trait,
            gene,
            genotype_values[:, column_index],
            outcome,
            covariate_values,
            continuous,
        )
        if result is not None:
            results.append(result)
    return results


def worker_count(requested: int | None, tasks: int) -> int:
    if tasks < 1:
        return 1
    if requested is None:
        return max(1, min(cpu_count() - 1, tasks))
    if requested < 1:
        raise ValueError("The number of CPU cores must be positive.")
    return requested


def run_observed_association(
    phenotype: pd.DataFrame,
    genotype: pd.DataFrame,
    covariate: pd.DataFrame,
    trait_types: dict[str, str],
    output: str,
    cores: int | None,
) -> None:
    phenotype, genotype, covariate = prepare_observed_inputs(phenotype, genotype, covariate)
    tasks = [
        (
            trait,
            phenotype[trait].to_numpy(),
            genotype,
            covariate.to_numpy(),
            trait_types,
        )
        for trait in phenotype.columns
    ]
    n_cores = worker_count(cores, len(tasks))
    start = time.time()
    print(f"Aligned samples: {len(phenotype)}")
    print(f"Covariates after preprocessing: {covariate.shape[1]}")
    print(
        f"Running observed association analysis for {len(tasks)} traits "
        f"with {n_cores} cores."
    )

    with Pool(n_cores) as pool:
        nested = list(
            progress(pool.imap(process_phenotype, tasks), len(tasks), "Observed traits")
        )

    rows = [row for group in nested for row in group]
    frame = pd.DataFrame(rows, columns=ASSOCIATION_COLUMNS)
    ensure_parent(output)
    frame.to_csv(output, sep=GENOTYPE_SEPARATOR, index=False)
    print(f"Saved {len(frame)} association rows to {output}.")
    print(f"Elapsed time: {time.time() - start:.1f} seconds.")


def worker_permutation_file(task: tuple[object, ...]) -> tuple[str, bool]:
    geno_path, output_path, phenotype, covariate_matrix, phenotype_index, trait_types = task
    try:
        genotype = read_matrix(str(geno_path), GENOTYPE_SEPARATOR).astype(float)
        common = phenotype_index.intersection(genotype.index)
        genotype = genotype.loc[common]
        phenotype = phenotype.loc[common]
        positions = phenotype_index.get_indexer(common)
        covariate_matrix = covariate_matrix[positions]
        genotype_values = genotype.to_numpy()
        rows = []

        for trait in phenotype.columns:
            continuous = trait_types[trait] == "quantitative"
            outcome, trait_genotypes, trait_covariates = phenotype_arrays(
                phenotype[trait].to_numpy(),
                genotype_values,
                covariate_matrix,
                continuous,
            )
            for column_index, gene in enumerate(genotype.columns):
                result = fit_pair(
                    trait,
                    gene,
                    trait_genotypes[:, column_index],
                    outcome,
                    trait_covariates,
                    continuous,
                )
                if result is not None and np.isfinite(result["P"]):
                    rows.append(
                        {
                            "Phenotype": trait,
                            "Genotype": gene,
                            "P": result["P"],
                        }
                    )

        ensure_parent(str(output_path))
        pd.DataFrame(rows, columns=["Phenotype", "Genotype", "P"]).to_csv(
            output_path, sep=GENOTYPE_SEPARATOR, index=False
        )
        return str(geno_path), True
    except Exception as error:
        print(f"Permutation failed for {geno_path}: {error}")
        return str(geno_path), False


def run_permutation_association(
    geno_dir: str,
    output_dir: str,
    phenotype: pd.DataFrame,
    covariate: pd.DataFrame,
    trait_types: dict[str, str],
    cores: int | None,
) -> None:
    genotype_files = sorted(
        [*Path(geno_dir).glob("*.csv"), *Path(geno_dir).glob("*.tsv")]
    )
    if not genotype_files:
        raise ValueError(f"No CSV or TSV genotype files found in {geno_dir}.")

    phenotype, covariate = prepare_permutation_inputs(phenotype, covariate)
    tasks = [
        (
            str(path),
            str(Path(output_dir) / f"res_{path.name}"),
            phenotype,
            covariate.to_numpy(),
            phenotype.index,
            trait_types,
        )
        for path in genotype_files
    ]
    n_cores = worker_count(cores, len(tasks))
    start = time.time()
    print(f"Found {len(tasks)} permutation genotype files.")
    print(
        f"Running association analysis for {len(tasks)} permutation datasets "
        f"with {n_cores} cores; this may take a while."
    )

    success = 0
    with Pool(n_cores) as pool:
        iterator = pool.imap(worker_permutation_file, tasks)
        for _, status in progress(iterator, len(tasks), "Permutation associations"):
            success += int(status)

    print(f"Completed {success}/{len(tasks)} permutation files.")
    print(f"Elapsed time: {(time.time() - start) / 60:.1f} minutes.")


def generate_permuted_genotypes(
    genotype_file: str,
    output_dir: str,
    n_perms: int = PERMUTATION_COUNT,
    seed: int = PERMUTATION_SEED,
) -> None:
    if n_perms < 1:
        raise ValueError("The number of permutations must be positive.")

    genotype = pd.read_csv(genotype_file, sep=GENOTYPE_SEPARATOR)
    if ID_COLUMN not in genotype.columns:
        raise ValueError(f"The genotype file must contain the '{ID_COLUMN}' column.")

    ensure_parent(Path(output_dir) / "placeholder")
    sample_ids = genotype[ID_COLUMN].to_numpy(copy=True)
    rng = np.random.default_rng(seed)
    print(f"Generating {n_perms} genotype permutations with seed {seed}.")

    for number in progress(
        range(1, n_perms + 1), n_perms, "Genotype permutations"
    ):
        permuted = genotype.copy()
        permuted[ID_COLUMN] = rng.permutation(sample_ids)
        output = Path(output_dir) / f"genotype_perm_{number:04d}.tsv"
        permuted.to_csv(output, sep=GENOTYPE_SEPARATOR, index=False)


def empirical_fdr(
    real_pvalues: np.ndarray,
    null_pvalues: list[np.ndarray],
) -> tuple[np.ndarray, float | None]:
    valid_real = ~np.isnan(real_pvalues)
    real = real_pvalues[valid_real]
    if len(real) == 0:
        return np.full(len(real_pvalues), np.nan), None

    thresholds = np.sort(np.unique(real))
    expected_null = np.zeros(len(thresholds), dtype=float)
    valid_permutations = 0
    for values in null_pvalues:
        values = values[~np.isnan(values)]
        if len(values) == 0:
            continue
        valid_permutations += 1
        expected_null += np.searchsorted(np.sort(values), thresholds, side="right")

    if valid_permutations == 0:
        raise ValueError("No valid permutation P values were found.")

    expected_null /= valid_permutations
    observed = np.searchsorted(np.sort(real), thresholds, side="right")
    estimated = np.minimum(expected_null / observed, 1.0)
    estimated = np.minimum.accumulate(estimated[::-1])[::-1]

    significant = np.flatnonzero(estimated <= TARGET_FDR)
    threshold = float(thresholds[significant[-1]]) if len(significant) else None
    adjusted = np.full(len(real_pvalues), np.nan)
    adjusted[valid_real] = estimated[np.searchsorted(thresholds, real, side="left")]
    print(f"Valid permutation runs used for eFDR: {valid_permutations}")
    return adjusted, threshold


def adjusted_threshold(pvalues: np.ndarray, adjusted: np.ndarray) -> float | None:
    significant = pvalues[adjusted <= TARGET_FDR]
    return float(np.max(significant)) if len(significant) else None


def threshold_text(value: float | None) -> str:
    return "NA" if value is None else f"{value:.2e}"


def exact_threshold_text(value: float | None) -> str:
    return "NA" if value is None else repr(float(value))


def run_multiple_testing(
    association_file: str,
    null_dir: str,
    output: str,
    trait_metadata: str = DEFAULT_TRAIT_METADATA,
    include_fdr_bh: bool = False,
) -> None:
    metadata = load_trait_metadata(trait_metadata)
    trait_types, _, _, _ = trait_metadata_maps(metadata)
    frame = pd.read_csv(association_file, sep=GENOTYPE_SEPARATOR)
    frame = frame.dropna(subset=["P"]).copy()
    if "N Carriers" in frame.columns:
        frame = frame[frame["N Carriers"] >= MIN_CARRIERS].copy()
    frame.reset_index(drop=True, inplace=True)

    validate_traits(frame["Trait"].dropna().unique(), trait_types)
    frame["Trait category"] = frame["Trait"].map(trait_types)

    pvalues = frame["P"].to_numpy(dtype=float)
    if include_fdr_bh:
        bh = stats.false_discovery_control(pvalues, method="bh")
        bh_threshold = adjusted_threshold(pvalues, bh)
        frame[f"FDR_BH ({threshold_text(bh_threshold)})"] = bh

    null_files = sorted([*Path(null_dir).glob("*.csv"), *Path(null_dir).glob("*.tsv")])
    if not null_files:
        raise ValueError(f"No permutation association files found in {null_dir}.")

    null_binary: list[np.ndarray] = []
    null_quantitative: list[np.ndarray] = []
    print(f"Loading null association results from {len(null_files)} files.")
    for path in progress(null_files, len(null_files), "Null result files"):
        null_frame = pd.read_csv(path, sep=GENOTYPE_SEPARATOR).dropna(subset=["P"])
        validate_traits(null_frame["Phenotype"].dropna().unique(), trait_types)
        categories = null_frame["Phenotype"].map(trait_types)
        null_binary.append(null_frame.loc[categories == "binary", "P"].to_numpy(dtype=float))
        null_quantitative.append(
            null_frame.loc[categories == "quantitative", "P"].to_numpy(dtype=float)
        )

    binary_mask = frame["Trait category"].eq("binary").to_numpy()
    quantitative_mask = frame["Trait category"].eq("quantitative").to_numpy()
    adjusted = np.full(len(frame), np.nan)
    binary_threshold = None
    quantitative_threshold = None

    if binary_mask.any():
        print("Calculating empirical FDR for binary traits.")
        adjusted_binary, binary_threshold = empirical_fdr(
            pvalues[binary_mask], null_binary
        )
        adjusted[binary_mask] = adjusted_binary
    if quantitative_mask.any():
        print("Calculating empirical FDR for quantitative traits.")
        adjusted_quantitative, quantitative_threshold = empirical_fdr(
            pvalues[quantitative_mask], null_quantitative
        )
        adjusted[quantitative_mask] = adjusted_quantitative

    print("Empirical FDR thresholds:")
    print(f"  Binary traits: {exact_threshold_text(binary_threshold)}")
    print(f"  Quantitative traits: {exact_threshold_text(quantitative_threshold)}")

    eFDR_column = (
        f"Empirical_FDR (Bin:{threshold_text(binary_threshold)}, "
        f"Cont:{threshold_text(quantitative_threshold)})"
    )
    frame[eFDR_column] = adjusted
    ensure_parent(output)
    frame.to_csv(output, sep=GENOTYPE_SEPARATOR, index=False)
    print(f"Saved multiple-testing results to {output}.")


def format_results(
    input_file: str,
    output: str,
    sort_by_p: bool,
    trait_metadata: str = DEFAULT_TRAIT_METADATA,
) -> None:
    frame = pd.read_csv(input_file, sep=GENOTYPE_SEPARATOR)
    metadata = load_trait_metadata(trait_metadata)
    trait_types, trait_labels, trait_domains, _ = trait_metadata_maps(metadata)
    validate_traits(frame["Trait"].dropna().unique(), trait_types)
    frame["Trait category"] = frame["Trait"].map(trait_types)
    eFDR_columns = [column for column in frame.columns if column.startswith("Empirical_FDR")]
    if not eFDR_columns:
        raise ValueError("No empirical FDR column found in the association results.")
    frame["Empirical_FDR_Sig"] = frame[eFDR_columns[0]].lt(TARGET_FDR).fillna(False)

    standard_columns = [
        "Trait",
        "Gene",
        "Beta",
        "SE",
        "OR",
        "95% CI Lower",
        "95% CI Upper",
        "P",
        "FDR_BH",
        "Empirical_FDR_Sig",
        "N Carriers",
        "Trait category",
        "Analysis method",
        "N Cases",
        "Affected cases",
        "Unaffected cases",
    ]
    column_map = {column.split(" (", 1)[0]: column for column in frame.columns}
    selected = [column_map[name] for name in standard_columns if name in column_map]
    frame = frame[selected].copy()
    frame.columns = [column.split(" (", 1)[0] for column in selected]
    frame["Clinical domain"] = frame["Trait"].map(trait_domains)
    frame["Trait"] = frame["Trait"].map(trait_labels)

    if sort_by_p:
        frame = frame.sort_values("P", ascending=True, na_position="last").reset_index(drop=True)
    ensure_parent(output)
    frame.to_csv(output, sep=GENOTYPE_SEPARATOR, index=False)
    print(f"Saved formatted results to {output}.")


def execute_association(args: argparse.Namespace) -> None:
    metadata = load_trait_metadata(args.trait_metadata)
    trait_types, _, _, _ = trait_metadata_maps(metadata)
    phenotype = read_matrix(args.pheno, PHENOTYPE_SEPARATOR)
    validate_traits(phenotype.columns, trait_types)
    covariate = read_matrix(args.cov, COVARIATE_SEPARATOR)
    if args.geno_dir:
        run_permutation_association(
            args.geno_dir,
            args.outdir,
            phenotype,
            covariate,
            trait_types,
            args.cores,
        )
        return
    if not args.geno:
        raise ValueError("--geno is required for observed association analysis.")
    genotype = read_matrix(args.geno, GENOTYPE_SEPARATOR)
    run_observed_association(
        phenotype,
        genotype,
        covariate,
        trait_types,
        args.output,
        args.cores,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="GPA reproducibility pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    association = subparsers.add_parser("association", help="Run observed or permutation association analysis")
    association.add_argument("--pheno", required=True)
    association.add_argument("--geno")
    association.add_argument("--cov", required=True)
    association.add_argument("--trait-metadata", default=DEFAULT_TRAIT_METADATA)
    association.add_argument("--output", default="geno_pheno_association_results.tsv")
    association.add_argument("--geno-dir")
    association.add_argument("--outdir", default="perm_results")
    association.add_argument("--cores", type=int, default=None)
    association.set_defaults(handler=execute_association)

    permutation = subparsers.add_parser("geno-permute", help="Generate genotype permutations")
    permutation.add_argument("--geno", required=True)
    permutation.add_argument("--outdir", required=True)
    permutation.add_argument("--n-perms", type=int, default=PERMUTATION_COUNT)
    permutation.add_argument("--seed", type=int, default=PERMUTATION_SEED)
    permutation.set_defaults(
        handler=lambda args: generate_permuted_genotypes(
            args.geno, args.outdir, args.n_perms, args.seed
        )
    )

    multitest = subparsers.add_parser(
        "multitest", help="Apply permutation eFDR correction"
    )
    multitest.add_argument("--assoc-file", required=True)
    multitest.add_argument("--null-pvals-dir", required=True)
    multitest.add_argument("--output", required=True)
    multitest.add_argument("--trait-metadata", default=DEFAULT_TRAIT_METADATA)
    multitest.add_argument(
        "--include-fdr-bh",
        action="store_true",
        help="Also include the Benjamini-Hochberg FDR column.",
    )
    multitest.set_defaults(
        handler=lambda args: run_multiple_testing(
            args.assoc_file,
            args.null_pvals_dir,
            args.output,
            args.trait_metadata,
            args.include_fdr_bh,
        )
    )

    formatter = subparsers.add_parser("format", help="Format corrected association results")
    formatter.add_argument("-i", "--input", required=True)
    formatter.add_argument("-o", "--output", required=True)
    formatter.add_argument("--trait-metadata", default=DEFAULT_TRAIT_METADATA)
    formatter.add_argument("--sort", action="store_true")
    formatter.set_defaults(
        handler=lambda args: format_results(
            args.input, args.output, args.sort, args.trait_metadata
        )
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.handler(args)


if __name__ == "__main__":
    main()
