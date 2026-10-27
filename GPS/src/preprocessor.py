import pickle

import numpy as np
import pandas as pd
from scipy.stats import norm

from .config import MISSINGNESS_THRESHOLD


class PhenotypePreprocessor:

    def __init__(
        self,
        continuous_traits,
        binary_traits,
        missing_threshold=MISSINGNESS_THRESHOLD,
    ):
        self.continuous_traits = list(continuous_traits)
        self.binary_traits = list(binary_traits)
        self.missing_threshold = float(missing_threshold)
        self.stats_mean = {}
        self.stats_prevalence = {}
        self.reference_values = {}
        self.fit_summary = {}
        self.is_fitted = False

    def filter_samples(self, df, missing_threshold=None):
        threshold = (
            self.missing_threshold
            if missing_threshold is None
            else float(missing_threshold)
        )
        traits = [
            trait for trait in self.continuous_traits + self.binary_traits
            if trait in df.columns
        ]
        if not traits:
            return df.copy()

        missing_rate = df[traits].isna().mean(axis=1)
        keep = missing_rate <= threshold
        removed = int((~keep).sum())
        print(
            f"Sample filter: removed {removed} samples with missingness "
            f"> {threshold:.1%}."
        )
        return df.loc[keep].copy()

    def _filter_samples(self, df):
        return self.filter_samples(df)

    @staticmethod
    def _rank_int(series):
        ranks = series.rank(na_option="keep")
        n = int(series.notna().sum())
        if n == 0:
            return pd.Series(np.nan, index=series.index, name=series.name)
        values = norm.ppf((ranks - 0.5) / n)
        return pd.Series(values, index=series.index, name=series.name)

    def _reference_rank_int(self, series, trait):
        reference = np.asarray(self.reference_values.get(trait, []), dtype=float)
        if reference.size == 0:
            raise ValueError(
                f"No derivation reference values are available for '{trait}'."
            )

        values = pd.to_numeric(series, errors="coerce")
        result = pd.Series(np.nan, index=series.index, dtype=float, name=series.name)
        valid = values.notna()
        if not valid.any():
            return result

        observed = values.loc[valid].to_numpy(dtype=float)
        left = np.searchsorted(reference, observed, side="left")
        right = np.searchsorted(reference, observed, side="right")
        average_rank = (left + right + 1.0) / 2.0
        n_reference = float(reference.size)
        probability = (average_rank - 0.5) / n_reference
        probability = np.clip(probability, 0.5 / n_reference, 1.0 - 0.5 / n_reference)
        result.loc[valid] = norm.ppf(probability)
        return result

    def fit_transform(self, df):
        clean = self.filter_samples(df)
        processed = clean.copy()

        for trait in self.continuous_traits:
            if trait in processed.columns:
                raw_values = pd.to_numeric(processed[trait], errors="coerce")
                self.reference_values[trait] = np.sort(
                    raw_values.dropna().to_numpy(dtype=float)
                )
                processed[trait] = self._rank_int(processed[trait])

        continuous_columns = [
            trait for trait in self.continuous_traits if trait in processed.columns
        ]
        binary_columns = [
            trait for trait in self.binary_traits if trait in processed.columns
        ]
        self.stats_mean = processed[continuous_columns].mean().to_dict()
        self.stats_prevalence = processed[binary_columns].mean().to_dict()
        self.fit_summary = {
            "input_n_samples": int(len(df)),
            "fit_n_samples": int(len(processed)),
            "removed_n_samples": int(len(df) - len(processed)),
            "missing_threshold": self.missing_threshold,
            "continuous_traits": self.continuous_traits,
            "binary_traits": self.binary_traits,
        }
        self.is_fitted = True
        return processed

    def transform(self, df, is_training=False, missing_threshold=None):
        if is_training:
            return self.fit_transform(df)
        if not self.is_fitted and not is_training:
            raise ValueError("The preprocessor must be fitted before transform().")

        clean = self.filter_samples(df, missing_threshold=missing_threshold)
        processed = clean.copy()

        for trait in self.continuous_traits:
            if trait in processed.columns:
                processed[trait] = self._reference_rank_int(processed[trait], trait)
        return processed

    def save(self, path):
        with open(path, "wb") as handle:
            pickle.dump(self, handle)
        print(f"Saved preprocessor: {path}")
