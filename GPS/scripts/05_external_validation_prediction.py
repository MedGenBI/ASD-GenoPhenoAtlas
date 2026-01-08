# scripts/05_external_validation_prediction.py

import sys
import os

import pandas as pd
import numpy as np
from itertools import product
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.config import (
    PROCESSED_DATA_DIR, RAW_DATA_PATH, OUTPUT_DIR, 
    CONTINUOUS_TRAITS, BINARY_TRAITS, 
    WGS_PHENO_PATH, WGS_PLP_PATH
)
from src.dataloader import load_training_data
from src.model import GPSModel

np.random.seed(42)

def find_best_global_params(df_pheno, df_assoc, df_plp, gene_list, preprocessor):

    print(f"Running Grid Search for {len(gene_list)} genes...")
    
    # Load data
    df_prev = pd.DataFrame(list(preprocessor.stats_prevalence.items()), columns=['phenotype', 'prevalence'])
    df_mean = pd.DataFrame(list(preprocessor.stats_mean.items()), columns=['phenotype', 'mean'])
    
    # Parameter grid
    pen_candidates = np.round(np.arange(0, 1.05, 0.05), 2)
    conf_candidates = np.round(np.arange(0, 1.05, 0.05), 2)
    grid = list(product(pen_candidates, conf_candidates))
    
    best_score = -1
    best_params = (0, 0)
    
    models_cache = []
    for gene in gene_list:
        gps_phes = df_assoc[(df_assoc['gene'] == gene) & (df_assoc['p'] < 0.05)]['phenotype'].tolist()
        if not gps_phes: continue

        model = GPSModel(gene, gps_phes, df_assoc, BINARY_TRAITS, CONTINUOUS_TRAITS)
        
        valid_mask = df_pheno[gps_phes].isna().mean(axis=1) <= 0.2
        if valid_mask.sum() == 0: continue
        
        df_valid = df_pheno[valid_mask].copy()

        carriers = df_plp[df_plp['Gene'] == gene]['spid'].unique()
        y_true = df_valid['spid'].isin(carriers).astype(int)
        
        if y_true.sum() == 0: continue
        
        # Cache models
        models_cache.append({
            'gene': gene,
            'model': model,
            'df': df_valid,
            'y_true': y_true
        })

    print(f"Prepared {len(models_cache)} valid gene models for grid search.")

    # Grid Search
    for idx, (pen, conf) in enumerate(grid):
        aucs = []
        for item in models_cache:
            model = item['model']
            df = item['df']
            y_true = item['y_true']

            model.fit(df, pen, conf, df_prev, df_mean)
            y_pred = model.predict(df, pen, conf, df_prev, df_mean)
            
            try:
                aucs.append(roc_auc_score(y_true, y_pred))
            except:
                pass
        
        if not aucs: continue
        mean_auc = np.mean(aucs)
        
        if mean_auc > best_score:
            best_score = mean_auc
            best_params = (pen, conf)
            
        if (idx + 1) % 20 == 0:
            print(f"Processed {idx+1}/{len(grid)} combos. Current best: {best_score:.4f} {best_params}")
            
    print(f"Best Params found: pen={best_params[0]}, conf={best_params[1]}, Mean AUC={best_score:.4f}")
    return best_params

def get_fitted_models(genes, df_pheno_deriv, df_assoc, pen, conf, preprocessor):
    """
    Get fitted GPSModel instances for the specified gene list

    :prams genes: list of genes
    :prams df_pheno_deriv: derivation cohort phenotype data (preprocessed)
    :prams df_assoc: genotype-phenotype association data
    :prams pen: penalty parameter
    :prams conf: confidence parameter
    :prams preprocessor: preprocessor object (used to retrieve statistics)

    :returns: dict, key=gene; value=fitted GPSModel
    """
    df_prev = pd.DataFrame(list(preprocessor.stats_prevalence.items()), columns=['phenotype', 'prevalence'])
    df_mean = pd.DataFrame(list(preprocessor.stats_mean.items()), columns=['phenotype', 'mean'])

    fitted_models = {}
    
    for gene in genes:
        gps_phes = df_assoc[(df_assoc['gene'] == gene) & (df_assoc['p'] < 0.05)]['phenotype'].tolist()
        if not gps_phes: continue
        
        model = GPSModel(gene, gps_phes, df_assoc, BINARY_TRAITS, CONTINUOUS_TRAITS)

        valid_mask = df_pheno_deriv[gps_phes].isna().mean(axis=1) <= 0.2
        if valid_mask.sum() == 0: continue
        df_valid = df_pheno_deriv[valid_mask].copy()

        model.fit(df_valid, pen, conf, df_prev, df_mean)
        
        fitted_models[gene] = model
    
    return fitted_models

def predict_wgs_cohort(df_wgs_phe, df_wgs_plp, gene_list, best_params, df_assoc, preprocessor, fitted_models):

    pen, conf = best_params
    df_prev = pd.DataFrame(list(preprocessor.stats_prevalence.items()), columns=['phenotype', 'prevalence'])
    df_mean = pd.DataFrame(list(preprocessor.stats_mean.items()), columns=['phenotype', 'mean'])
    
    results = []
    
    for gene in gene_list:

        gene_plp = df_wgs_plp[df_wgs_plp['Gene'] == gene]
        if gene_plp.empty: continue
        
        carriers_spids = gene_plp['spid'].unique()

        gps_phes = df_assoc[(df_assoc['gene'] == gene) & (df_assoc['p'] < 0.05)]['phenotype'].tolist()
        if not gps_phes: continue
    
        # Get the fitted model
        model = fitted_models.get(gene)
        
        # Dynamic data filtering (missing rate <= 20%)
        valid_mask = df_wgs_phe[gps_phes].isna().mean(axis=1) <= 0.2
        df_valid = df_wgs_phe[valid_mask].copy()
        
        y_true = df_valid['spid'].isin(carriers_spids).astype(int)
        if y_true.sum() == 0: continue

        scores = model.predict(df_valid, pen, conf, df_prev, df_mean)
        
        res_df = pd.DataFrame({
            'spid': df_valid['spid'],
            f'{gene}_score': scores,
            f'{gene}_carrier': y_true
        })
        results.append(res_df.set_index('spid'))
        
    if not results:
        return pd.DataFrame()
        
    return pd.concat(results, axis=1)

def main():
    print("External Validation - Prediction...")
    
    # 1. Load derivation cohort data
    df_pheno_deriv, _, df_assoc, df_plp_deriv, preprocessor = load_training_data(RAW_DATA_PATH, PROCESSED_DATA_DIR)
    
    # 2. Load and preprocess WGS data
    print("\nLoading WGS data...")
    df_wgs_phe = pd.read_csv(WGS_PHENO_PATH)
    df_wgs_plp = pd.read_csv(WGS_PLP_PATH, sep='\t')
    
    print("\nPreprocessing WGS data (INT using Derivation ECDF)...")
    
    # 2.1 Missing rate filtering
    missing_rate = 0.2
    cols_to_check = CONTINUOUS_TRAITS + BINARY_TRAITS
    valid_cols = [c for c in cols_to_check if c in df_wgs_phe.columns]
    df_wgs_phe = df_wgs_phe[df_wgs_phe[valid_cols].isna().mean(axis=1) <= missing_rate].copy()
    
    # 2.2 INT transformation
    # Get ECDF from preprocessor object
    for pheno in CONTINUOUS_TRAITS:
        if pheno in df_wgs_phe.columns and pheno in preprocessor.ecdfs:
            ecdf = preprocessor.ecdfs[pheno]
            original_values = df_wgs_phe[pheno].dropna()

            if len(original_values) > 0:
                quantiles = ecdf(original_values)
                quantiles = np.clip(quantiles, 1e-6, 1 - 1e-6)

                rint_values = norm.ppf(quantiles)
                df_wgs_phe.loc[original_values.index, pheno] = rint_values
                
    print(f"WGS Phenotypes prepared. Shape: {df_wgs_phe.shape}")
    
    # 3. Grouping
    gene_metrics = pd.read_csv(os.path.join(PROCESSED_DATA_DIR, 'gene_level_evaluation_metrics.tsv'), sep='\t')
    hc_genes = gene_metrics[(gene_metrics['auc_median'] >= 0.8) & (gene_metrics['fdr_phenotypes'] > 0)]['gene'].tolist()
    ex_genes = gene_metrics[gene_metrics['auc_median'] >= 0.8]['gene'].tolist()
    
    # 4. High-confidence Set
    print("\n--- Processing High-confidence Set ---")

    best_params_hc = find_best_global_params(df_pheno_deriv, df_assoc, df_plp_deriv, hc_genes, preprocessor)
    fitted_models = get_fitted_models(hc_genes, df_pheno_deriv, df_assoc, best_params_hc[0], best_params_hc[1], preprocessor)
    df_pred_wgs_hc = predict_wgs_cohort(df_wgs_phe, df_wgs_plp, hc_genes, best_params_hc, df_assoc, preprocessor, fitted_models=fitted_models)
    
    # 5. Extended Set
    print("\n--- Processing Extended Set ---")

    best_params_ex = find_best_global_params(df_pheno_deriv, df_assoc, df_plp_deriv, ex_genes, preprocessor)
    fitted_models = get_fitted_models(ex_genes, df_pheno_deriv, df_assoc, best_params_ex[0], best_params_ex[1], preprocessor)
    df_pred_wgs_ex = predict_wgs_cohort(df_wgs_phe, df_wgs_plp, ex_genes, best_params_ex, df_assoc, preprocessor, fitted_models=fitted_models)
    
    # 6. Save predictions
    out_hc = os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_predictions_hc.tsv')
    out_ex = os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_predictions_ex.tsv')
    
    df_pred_wgs_hc.to_csv(out_hc, sep='\t')
    df_pred_wgs_ex.to_csv(out_ex, sep='\t')
    
    print("WGS Prediction completed.")


if __name__ == "__main__":
    main()
