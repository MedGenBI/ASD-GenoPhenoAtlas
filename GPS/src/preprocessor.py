# src/preprocessor.py

import pickle
from typing import List, Dict

import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF


class PhenotypePreprocessor:
    def __init__(self, continuous_traits: List[str], binary_traits: List[str], missing_threshold: float = 0.2):

        self.continuous_traits = continuous_traits
        self.binary_traits = binary_traits
        self.missing_threshold = missing_threshold

        self.ecdfs: Dict[str, ECDF] = {}                # Save ECDF
        self.stats_mean: Dict[str, float] = {}          # Save mean of continuous traits
        self.stats_prevalence: Dict[str, float] = {}    # Save prevalence of binary traits
        self.is_fitted = False

    @staticmethod
    def _rank_based_int(series: pd.Series) -> pd.Series:

        ranks = series.rank(na_option='keep')
        n = series.notna().sum()
        transformed = norm.ppf((ranks - 0.5) / n)

        return pd.Series(transformed, index=series.index, name=series.name)

    def _filter_samples(self, df: pd.DataFrame) -> pd.DataFrame:
        
        target_cols = self.continuous_traits + self.binary_traits
        valid_cols = [c for c in target_cols if c in df.columns]
        
        if not valid_cols:
            return df
            
        missing_rates = df[valid_cols].isna().mean(axis=1)
        kept_mask = missing_rates <= self.missing_threshold
        
        n_removed = (~kept_mask).sum()
        print(f"Filtering: Removed {n_removed} samples with missing rate > {self.missing_threshold:.1%}")
        
        return df[kept_mask].copy()

    def fit(self, df: pd.DataFrame):

        print("Fitting preprocessor: Learning ECDFs...")

        df_clean = self._filter_samples(df)
        
        for trait in self.continuous_traits:
            if trait in df_clean.columns:
                valid_values = df_clean[trait].dropna()
                if len(valid_values) > 0:
                    self.ecdfs[trait] = ECDF(valid_values)
        
        self.is_fitted = True
        return self

    def transform(self, df: pd.DataFrame, is_training: bool = True) -> pd.DataFrame:
        """
        Data transform
        """
        if not self.is_fitted:
            raise ValueError("Preprocessor must be fitted before transform.")

        # 1. Samples filtering
        df_processed = self._filter_samples(df).copy()
        
        # 2. INT transformation
        print("Applying Rank-based Inverse Normal Transformation...")
        for trait in self.continuous_traits:
            if trait in df_processed.columns:
                df_processed[trait] = self._rank_based_int(df_processed[trait])

        if is_training:
            print("Calculating statistics (Mean & Prevalence) on transformed data...")
            # Continuous traits mean
            self.stats_mean = df_processed[self.continuous_traits].mean().to_dict()
            # Binary traits prevalence
            self.stats_prevalence = df_processed[self.binary_traits].mean().to_dict()
            
        return df_processed

    def fit_transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """fit & transform"""
        self.fit(df)
        return self.transform(df, is_training=True)
    
    def save(self, path: str):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        print(f"Preprocessor saved to {path}")
