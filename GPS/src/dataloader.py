import os
import pickle

import pandas as pd

from .config import (
    ASSOCIATION_P_THRESHOLD,
    FDR_SIG_COL,
    GENE_COL,
    P_COL,
    PHENO_COL,
    PHENO_SAMPLE_COL,
    normalize_trait_id,
)
from .preprocessor import PhenotypePreprocessor


GENE_ASSOCIATION_COUNTS_FILENAME = "gene_association_counts.csv"
PHENOTYPE_ASSOCIATION_COUNTS_FILENAME = "phenotype_association_counts.csv"
PROCESSED_PHENOTYPE_FILENAME = "phenotypes_processed_int.csv"
PREPROCESSOR_FILENAME = "preprocessor_model.pkl"


def read_table(path, index_col=None):
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    lower = path.lower()
    if lower.endswith(".csv"):
        return pd.read_csv(path, index_col=index_col)
    if lower.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", index_col=index_col)
    raise ValueError(f"Unsupported input format: {path}")


def normalize_association_traits(association):
    """Map association Trait values to the canonical trait_id namespace."""
    if PHENO_COL not in association.columns:
        raise ValueError(f"Association table must contain {PHENO_COL}.")

    result = association.copy()
    unknown = set()
    normalized = []
    for value in result[PHENO_COL]:
        if pd.isna(value) or not str(value).strip():
            normalized.append(value)
            continue
        try:
            normalized.append(normalize_trait_id(value))
        except ValueError:
            unknown.add(str(value).strip())
            normalized.append(value)
    if unknown:
        values = ", ".join(sorted(unknown)[:10])
        suffix = "..." if len(unknown) > 10 else ""
        raise ValueError(
            "Association table contains trait names not defined in trait metadata: "
            f"{values}{suffix}"
        )
    result[PHENO_COL] = normalized
    return result


def normalize_phenotype_columns(phenotype):
    """Rename phenotype columns to trait_id while preserving the sample column."""
    if PHENO_SAMPLE_COL not in phenotype.columns:
        raise ValueError(
            f"Phenotype table must contain the sample column {PHENO_SAMPLE_COL}."
        )

    rename = {}
    unknown = []
    for column in phenotype.columns:
        if column == PHENO_SAMPLE_COL:
            continue
        try:
            canonical = normalize_trait_id(column)
        except ValueError:
            unknown.append(str(column))
            continue
        rename[column] = canonical

    if unknown:
        preview = ", ".join(unknown[:10])
        suffix = "..." if len(unknown) > 10 else ""
        raise ValueError(
            "Phenotype table contains trait columns not defined in trait metadata: "
            f"{preview}{suffix}"
        )

    result = phenotype.rename(columns=rename).copy()
    duplicate_columns = result.columns[result.columns.duplicated()].tolist()
    if duplicate_columns:
        duplicates = ", ".join(dict.fromkeys(map(str, duplicate_columns)))
        raise ValueError(
            "Phenotype table contains duplicate trait columns after trait-name "
            f"normalization: {duplicates}"
        )
    return result


def load_raw_data(association_path, phenotype_path):
    print(f"Loading association data: {association_path}")
    association = normalize_association_traits(read_table(association_path))
    print(f"Loading phenotype data: {phenotype_path}")
    phenotype = normalize_phenotype_columns(read_table(phenotype_path))
    return None, association, phenotype


def generate_association_count_summaries(
    association,
    p_col=P_COL,
    fdr_sig_col=FDR_SIG_COL,
    gene_col=GENE_COL,
    pheno_col=PHENO_COL,
):
    if not {p_col, fdr_sig_col, gene_col, pheno_col}.issubset(association.columns):
        missing = sorted(
            {p_col, fdr_sig_col, gene_col, pheno_col} - set(association.columns)
        )
        raise ValueError(f"Association table is missing columns: {missing}")

    p_values = pd.to_numeric(association[p_col], errors="coerce")
    fdr_values = association[fdr_sig_col].astype(bool)
    work = association.copy()
    work["_nominal"] = p_values < ASSOCIATION_P_THRESHOLD
    work["_fdr"] = fdr_values

    gene_counts = work.groupby(gene_col).agg(
        p_phenotypes=("_nominal", "sum"),
        fdr_phenotypes=("_fdr", "sum"),
    )
    gene_counts = gene_counts.sort_values(
        ["fdr_phenotypes", "p_phenotypes"],
        ascending=False,
    )

    phenotype_counts = work.groupby(pheno_col).agg(
        p_genes=("_nominal", "sum"),
        fdr_genes=("_fdr", "sum"),
    )
    phenotype_counts = phenotype_counts.sort_values(
        "fdr_genes", ascending=False,
    )
    return gene_counts, phenotype_counts


def load_training_data(assoc_path, plp_path, processed_data_dir):
    phenotype_path = os.path.join(
        processed_data_dir, PROCESSED_PHENOTYPE_FILENAME
    )
    gene_counts_path = os.path.join(
        processed_data_dir, GENE_ASSOCIATION_COUNTS_FILENAME
    )
    preprocessor_path = os.path.join(
        processed_data_dir, PREPROCESSOR_FILENAME
    )

    phenotype = normalize_phenotype_columns(read_table(phenotype_path))
    gene_counts = read_table(gene_counts_path, index_col=GENE_COL)
    association = normalize_association_traits(read_table(assoc_path))
    plp = read_table(plp_path)

    with open(preprocessor_path, "rb") as handle:
        preprocessor = pickle.load(handle)

    print("Loaded derivation data.")
    return phenotype, gene_counts, association, plp, preprocessor
