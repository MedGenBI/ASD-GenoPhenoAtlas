import numpy as np
import pandas as pd

from .config import BETA_COL, FDR_SIG_COL, GENE_COL, PHENO_COL


class GPSModel:

    def __init__(
        self,
        target_gene,
        gene_assoc_phenos,
        association,
        binary_traits,
        continuous_traits,
        continuous_traits_direction,
    ):
        self.gene = str(target_gene)
        self.gene_assoc_phenos = list(gene_assoc_phenos)
        self.binary_traits = set(binary_traits)
        self.continuous_traits = set(continuous_traits)
        self.direction = {
            str(trait): str(value).strip().lower()
            for trait, value in continuous_traits_direction.items()
        }
        self.effects = self._extract_effects(association)
        self.normalization_params = {}
        self.is_fitted = False

    def _extract_effects(self, association):
        rows = association.loc[
            (association[GENE_COL].astype(str) == self.gene)
            & association[PHENO_COL].isin(self.gene_assoc_phenos)
        ]
        effects = {}
        for _, row in rows.iterrows():
            beta = pd.to_numeric(pd.Series([row[BETA_COL]]), errors="coerce").iloc[0]
            if not np.isfinite(beta):
                raise ValueError(f"Invalid beta for {self.gene}:{row[PHENO_COL]}")
            effects[str(row[PHENO_COL])] = {
                "beta": float(beta),
                "fdr": bool(row[FDR_SIG_COL]),
            }
        return effects

    def _abnormal_deviation(self, trait, values, mean_value):
        direction = self.direction.get(trait)
        if direction in {"lower", "decrease", "low"}:
            return np.maximum(0.0, mean_value - values)
        if direction in {"higher", "later", "increase", "high"}:
            return np.maximum(0.0, values - mean_value)
        raise ValueError(
            f"Unsupported direction '{direction}' for continuous trait '{trait}'."
        )

    def _raw_components(self, phenotype, penalty, confidence, prevalence, means):
        score_binary = pd.Series(0.0, index=phenotype.index)
        score_continuous = pd.Series(0.0, index=phenotype.index)

        for trait in self.gene_assoc_phenos:
            if trait not in phenotype.columns or trait not in self.effects:
                continue
            effect = self.effects[trait]
            values = phenotype[trait].fillna(0)
            observed = phenotype[trait].notna().astype(float)
            weight = 1.0 if effect["fdr"] else confidence
            contribution = values * effect["beta"] * weight * observed

            if trait in self.binary_traits:
                score_binary += contribution
            elif trait in self.continuous_traits:
                score_continuous += contribution

        prevalence_map = prevalence.set_index("phenotype")["prevalence"].to_dict()
        for trait in self.binary_traits.difference(self.gene_assoc_phenos):
            if trait not in phenotype.columns or trait not in prevalence_map:
                continue
            observed = (phenotype[trait] == 1) & phenotype[trait].notna()
            score_binary -= (
                (1.0 - float(prevalence_map[trait]))
                * observed.astype(float)
                * penalty
            )

        mean_map = means.set_index("phenotype")["mean"].to_dict()
        for trait in self.continuous_traits.difference(self.gene_assoc_phenos):
            if trait not in phenotype.columns or trait not in mean_map:
                continue
            values = phenotype[trait]
            deviation = self._abnormal_deviation(
                trait, values, float(mean_map[trait])
            )
            score_continuous -= (
                np.where(values.notna(), deviation, 0.0) * penalty
            )

        return score_binary, score_continuous

    def fit(self, phenotype, penalty, confidence, prevalence, means):
        score_binary, score_continuous = self._raw_components(
            phenotype, penalty, confidence, prevalence, means
        )
        binary_std = score_binary.std()
        continuous_std = score_continuous.std()
        self.normalization_params = {
            "mean_binary": float(score_binary.mean()),
            "std_binary": float(binary_std) if binary_std > 1e-9 else 1.0,
            "mean_continuous": float(score_continuous.mean()),
            "std_continuous": (
                float(continuous_std) if continuous_std > 1e-9 else 1.0
            ),
        }
        self.is_fitted = True
        return self

    def predict(self, phenotype, penalty, confidence, prevalence, means):
        if not self.is_fitted:
            raise ValueError("Call fit before predict.")

        score_binary, score_continuous = self._raw_components(
            phenotype, penalty, confidence, prevalence, means
        )
        params = self.normalization_params
        z_binary = (
            score_binary - params["mean_binary"]
        ) / params["std_binary"]
        z_continuous = (
            score_continuous - params["mean_continuous"]
        ) / params["std_continuous"]
        final = (z_binary + z_continuous) / 2.0
        result = 1.0 / (1.0 + np.exp(-final))
        return result.fillna(0.0)
