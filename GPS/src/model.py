# src/model.py

import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Tuple

class GPSModel:
    def __init__(
        self, 
        target_gene: str,
        gps_phenotypes: List[str],
        assoc_df: pd.DataFrame,
        binary_traits: List[str],
        continuous_traits: List[str],
        k: float = 1.0,
        x0: float = 0.0
    ):
        """
        Initialize GPS model
        
        :param target_gene: target gene name
        :param gps_phenotypes: list of phenotypes nominally associated with the gene (used for scoring)
        :param assoc_df: association results DataFrame (containing gene, phenotype, beta, p_fdr columns)
        :param binary_traits: list of all binary phenotype names
        :param continuous_traits: list of all continuous phenotype names
        """
        self.gene = target_gene
        self.gps_phenotypes = gps_phenotypes
        self.binary_traits = set(binary_traits)
        self.continuous_traits = set(continuous_traits)
        self.k = k
        self.x0 = x0
        
        self.effects = self._extract_effects(assoc_df)
        
        self.normalization_params = {} 
        self.is_fitted = False

    def _extract_effects(self, assoc_df: pd.DataFrame) -> Dict:
        gene_assoc = assoc_df[assoc_df['gene'] == self.gene]
        effects = {}
        for _, row in gene_assoc.iterrows():
            if row['phenotype'] in self.gps_phenotypes:
                effects[row['phenotype']] = {
                    'beta': row.get('beta', 0),
                    'p_fdr': row.get('p_fdr', 1.0)
                }
        return effects

    def _compute_raw_scores(
        self, 
        df_pheno: pd.DataFrame, 
        pen_factor: float, 
        conf_factor: float,
        df_prev: Optional[pd.DataFrame] = None,
        df_mean: Optional[pd.DataFrame] = None
    ) -> Tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
        """
        Compute raw scores (un-normalized)
        """
        score_bin = pd.Series(0.0, index=df_pheno.index)
        score_cont = pd.Series(0.0, index=df_pheno.index)
        
        # 1. Compute contributions
        for phe in self.gps_phenotypes:
            if phe not in df_pheno.columns: continue
            if phe not in self.effects: continue
            
            effect_data = self.effects[phe]
            beta = effect_data['beta']
            weight = 1.0 if effect_data['p_fdr'] < 0.05 else conf_factor
            
            vals = df_pheno[phe].fillna(0)
            mask = df_pheno[phe].notna().astype(float)
            
            contribution = vals * beta * weight * mask
            
            if phe in self.binary_traits:
                score_bin += contribution
            elif phe in self.continuous_traits:
                score_cont += contribution
        
        # 2. Compute penalty terms - binary traits
        if df_prev is not None and pen_factor > 0:
            # find traits not in gps_phenotypes
            penalty_traits = [p for p in self.binary_traits if p not in self.gps_phenotypes]
            prev_dict = df_prev.set_index('phenotype')['prevalence'].to_dict()
            
            for phe in penalty_traits:
                if phe in df_pheno.columns and phe in prev_dict:
                    prev = prev_dict[phe]
                    mask = (df_pheno[phe] == 1) & df_pheno[phe].notna()
                    score_bin -= (1 - prev) * mask.astype(float) * pen_factor

        # 3. Compute penalty terms - continuous traits
        if df_mean is not None and pen_factor > 0:
            penalty_traits = [p for p in self.continuous_traits if p not in self.gps_phenotypes]
            mean_dict = df_mean.set_index('phenotype')['mean'].to_dict()
            
            for phe in penalty_traits:
                if phe in df_pheno.columns and phe in mean_dict:
                    mean_val = mean_dict[phe]
                    vals = df_pheno[phe]
                    mask = vals.notna()

                    if phe.startswith(('DCDQ', 'ABC')):
                         diff = np.maximum(0, mean_val - vals)
                    else:
                         diff = np.maximum(0, vals - mean_val)
                    score_cont -= np.where(mask, diff, 0) * pen_factor

        return score_bin, score_cont

    def fit(
        self, 
        df_pheno: pd.DataFrame, 
        pen_factor: float, 
        conf_factor: float,
        df_prev: pd.DataFrame,
        df_mean: pd.DataFrame
    ):
        s_bin, s_cont = self._compute_raw_scores(df_pheno, pen_factor, conf_factor, df_prev, df_mean)

        self.normalization_params = {
            'mean_bin': s_bin.mean(),
            'std_bin': s_bin.std() if s_bin.std() > 1e-9 else 1.0,
            'mean_cont': s_cont.mean(),
            'std_cont': s_cont.std() if s_cont.std() > 1e-9 else 1.0
        }
        self.is_fitted = True
        return self

    def predict(
        self, 
        df_pheno: pd.DataFrame, 
        pen_factor: float, 
        conf_factor: float,
        df_prev: pd.DataFrame,
        df_mean: pd.DataFrame
    ) -> pd.Series:

        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")
            
        s_bin, s_cont = self._compute_raw_scores(df_pheno, pen_factor, conf_factor, df_prev, df_mean)
        
        params = self.normalization_params
        
        # Z-score standardization
        z_bin = (s_bin - params['mean_bin']) / params['std_bin']
        z_cont = (s_cont - params['mean_cont']) / params['std_cont']
        
        final_score = (z_bin + z_cont) / 2.0
        
        # Sigmoid
        final_prob = 1.0 / (1.0 + np.exp(-self.k * final_score + self.k * self.x0))
        
        return final_prob.fillna(0)
    