import csv
import os

_CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.abspath(os.path.join(_CONFIG_DIR, "../../Data"))
DATA_DIR = os.environ.get("GPS_DATA_DIR", DATA_DIR)

PHENOTYPE_PATH = os.environ.get(
    "GPS_PHENOTYPE_PATH", os.path.join(DATA_DIR, "phenotype.csv")
)
PLP_PATH = os.environ.get(
    "GPS_PLP_PATH", os.path.join(DATA_DIR, "ACMG_PLP_variants.tsv")
)
OUTPUT_DIR = os.environ.get("GPS_OUTPUT_DIR", "./output")
WGS_PHENO_PATH = os.environ.get(
    "GPS_WGS_PHENO_PATH", os.path.join(DATA_DIR, "phenotype_GS.csv")
)
WGS_PLP_PATH = os.environ.get(
    "GPS_WGS_PLP_PATH", os.path.join(DATA_DIR, "ACMG_PLP_variants_GS.tsv")
)
TRAIT_METADATA = os.environ.get(
    "GPS_TRAIT_METADATA", os.path.join(DATA_DIR, "trait_metadata.csv")
)

P_COL = "P"
BETA_COL = "Beta"
FDR_SIG_COL = "Empirical_FDR_Sig"
GENE_COL = "Gene"
PHENO_COL = "Trait"
PLP_SAMPLE_COL = "Sample ID"
PHENO_SAMPLE_COL = "spid"

ASSOCIATION_P_THRESHOLD = 0.05
MISSINGNESS_THRESHOLD = 0.20


def _read_trait_metadata(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Trait metadata file not found: {path}. "
            "Set GPS_TRAIT_METADATA to the metadata CSV."
        )

    required = {
        "trait_id",
        "trait_label",
        "trait_type",
        "clinical_domain",
        "effect_direction",
    }
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Trait metadata file has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                "Trait metadata is missing required columns: "
                + ", ".join(sorted(missing))
            )

        records = []
        seen = set()
        for row_number, row in enumerate(reader, start=2):
            trait_id = str(row["trait_id"] or "").strip()
            trait_label = str(row["trait_label"] or "").strip()
            trait_type = str(row["trait_type"] or "").strip().lower()
            clinical_domain = str(row["clinical_domain"] or "").strip()
            direction = str(row["effect_direction"] or "").strip().lower()

            if not trait_id:
                raise ValueError(f"Empty trait_id at metadata row {row_number}.")
            if not trait_label:
                raise ValueError(f"Empty trait_label at metadata row {row_number}.")
            if trait_id in seen:
                raise ValueError(f"Duplicate trait_id in metadata: {trait_id}")
            if trait_type not in {"binary", "quantitative", "continuous"}:
                raise ValueError(
                    f"Unsupported trait_type '{trait_type}' for {trait_id}."
                )
            if not clinical_domain:
                raise ValueError(f"Empty clinical_domain for {trait_id}.")

            if trait_type in {"quantitative", "continuous"}:
                if direction not in {"higher", "later", "lower"}:
                    raise ValueError(
                        f"Unsupported effect_direction '{direction}' for {trait_id}."
                    )
                # Developmental milestones use "later", which is higher severity.
                direction = "higher" if direction == "later" else direction
            else:
                direction = ""

            records.append(
                {
                    "trait_id": trait_id,
                    "trait_label": trait_label,
                    "trait_type": (
                        "quantitative" if trait_type == "continuous" else trait_type
                    ),
                    "clinical_domain": clinical_domain,
                    "effect_direction": direction,
                }
            )
            seen.add(trait_id)

    if not records:
        raise ValueError(f"Trait metadata file is empty: {path}")
    return records


_TRAIT_RECORDS = _read_trait_metadata(TRAIT_METADATA)
TRAIT_ID_TO_LABEL = {row["trait_id"]: row["trait_label"] for row in _TRAIT_RECORDS}
TRAIT_NAME_TO_ID = {}
for _record in _TRAIT_RECORDS:
    _trait_id = _record["trait_id"]
    for _trait_name in (_trait_id, _record["trait_label"]):
        _existing_id = TRAIT_NAME_TO_ID.get(_trait_name)
        if _existing_id is not None and _existing_id != _trait_id:
            raise ValueError(
                "Trait metadata contains an ambiguous trait name: "
                f"{_trait_name} maps to both {_existing_id} and {_trait_id}."
            )
        TRAIT_NAME_TO_ID[_trait_name] = _trait_id
TRAIT_TYPES = {row["trait_id"]: row["trait_type"] for row in _TRAIT_RECORDS}
TRAIT_DOMAINS = {
    row["trait_id"]: row["clinical_domain"] for row in _TRAIT_RECORDS
}


def normalize_trait_id(value):
    """Return the canonical trait_id for an ID or display label."""
    name = str(value).strip()
    if not name:
        raise ValueError("Trait name cannot be empty.")
    try:
        return TRAIT_NAME_TO_ID[name]
    except KeyError as error:
        raise ValueError(
            f"Unknown trait name '{name}'. It is not defined in {TRAIT_METADATA}."
        ) from error


def _traits_in_domain(domain):
    return [
        row["trait_id"] for row in _TRAIT_RECORDS
        if row["clinical_domain"] == domain
    ]


SCALE = _traits_in_domain("Scales")
DEV_MILESTONES = _traits_in_domain("Developmental Milestones")
MH_PC = _traits_in_domain("Mental Health and Psychiatric Conditions")
NDD_RELATED = _traits_in_domain("Neurodevelopmental and Related Conditions")
PERINATAL = _traits_in_domain("Perinatal Complications")
BIRTH_DEFECTS = _traits_in_domain("Birth Defects")
GROWTH = _traits_in_domain("Growth Abnormalities")

# Preserve the historical scale subgroups while deriving them from metadata.
SCQ = [trait for trait in SCALE if trait.startswith("SCQ_")]
RBSR = [trait for trait in SCALE if trait.startswith("RBSR_")]
DCDQ = [trait for trait in SCALE if trait.startswith("DCDQ_")]
ABC = [trait for trait in SCALE if trait.startswith("ABC_")]

CONTINUOUS_TRAITS = [
    row["trait_id"] for row in _TRAIT_RECORDS
    if row["trait_type"] == "quantitative"
]
BINARY_TRAITS = [
    row["trait_id"] for row in _TRAIT_RECORDS
    if row["trait_type"] == "binary"
]
ALL_TRAITS = [row["trait_id"] for row in _TRAIT_RECORDS]
CONTINUOUS_TRAITS_DIRECTION = {
    row["trait_id"]: row["effect_direction"]
    for row in _TRAIT_RECORDS
    if row["trait_type"] == "quantitative"
}
