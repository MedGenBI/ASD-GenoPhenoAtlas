# scripts/02_model_optimization.py

import sys
import os
import pandas as pd
import numpy as np
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.metrics import roc_auc_score
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.config import (
    PROCESSED_DATA_DIR, RAW_DATA_PATH, 
    BINARY_TRAITS, CONTINUOUS_TRAITS
)
from src.model import GPSModel
from src.dataloader import load_training_data


GLOBAL_DATA = {
    'df_pheno': None,
    'df_assoc': None,
    'df_plp': None,
    'df_prev': None,
    'df_mean': None
}

def init_worker(pheno, assoc, plp, prev, mean):

    GLOBAL_DATA['df_pheno'] = pheno
    GLOBAL_DATA['df_assoc'] = assoc
    GLOBAL_DATA['df_plp'] = plp
    GLOBAL_DATA['df_prev'] = prev
    GLOBAL_DATA['df_mean'] = mean

def process_one_gene_log_cv(test_gene, all_genes, param_grid, df_pheno=None, df_assoc=None, df_plp=None, df_prev=None, df_mean=None):

    df_pheno = GLOBAL_DATA['df_pheno'] if df_pheno is None else df_pheno
    df_assoc = GLOBAL_DATA['df_assoc'] if df_assoc is None else df_assoc
    df_plp = GLOBAL_DATA['df_plp'] if df_plp is None else df_plp
    df_prev = GLOBAL_DATA['df_prev'] if df_prev is None else df_prev
    df_mean = GLOBAL_DATA['df_mean'] if df_mean is None else df_mean

    # 1. Training genes
    train_genes = [g for g in all_genes if g != test_gene]
    
    best_inner_score = -1.0
    best_params = (0, 0)

    # Grid Search on Train Genes
    for pen, conf in param_grid:
        fold_aurocs = []

        for train_gene in train_genes:
            # 1. Prepare model
            gps_phes = df_assoc[(df_assoc['gene'] == train_gene) & (df_assoc['p'] < 0.05)]['phenotype'].tolist()
            if not gps_phes: continue
            
            model = GPSModel(train_gene, gps_phes, df_assoc, BINARY_TRAITS, CONTINUOUS_TRAITS)
            
            # 2. Prepare data
            valid_mask = df_pheno[gps_phes].isna().mean(axis=1) <= 0.2
            df_valid = df_pheno[valid_mask].copy()
            
            carriers = df_plp[df_plp['Gene'] == train_gene]['spid'].unique()
            y_true = df_valid['spid'].isin(carriers).astype(int)
            
            if y_true.sum() == 0: continue
            
            # 3. Fit & Predict
            model.fit(df_valid, pen, conf, df_prev, df_mean)
            y_pred = model.predict(df_valid, pen, conf, df_prev, df_mean)
            
            try:
                auc_val = roc_auc_score(y_true, y_pred)
                fold_aurocs.append(auc_val)
            except:
                pass
        
        mean_auc = np.mean(fold_aurocs) if fold_aurocs else 0
        
        if mean_auc > best_inner_score:
            best_inner_score = mean_auc
            best_params = (pen, conf)

    # Predict on Test Gene
    pen_opt, conf_opt = best_params
    
    gps_phes_test = df_assoc[(df_assoc['gene'] == test_gene) & (df_assoc['p'] < 0.05)]['phenotype'].tolist()
    
    result_data = {
        'gene': test_gene,
        'best_pen': pen_opt,
        'best_conf': conf_opt,
        'cv_train_auc': best_inner_score,
        'pred_df': None
    }
    
    if gps_phes_test:

        model_test = GPSModel(test_gene, gps_phes_test, df_assoc, BINARY_TRAITS, CONTINUOUS_TRAITS)

        valid_mask = df_pheno[gps_phes_test].isna().mean(axis=1) <= 0.2
        df_valid_test = df_pheno[valid_mask].copy()

        model_test.fit(df_valid_test, pen_opt, conf_opt, df_prev, df_mean)

        y_pred_test = model_test.predict(df_valid_test, pen_opt, conf_opt, df_prev, df_mean)

        carriers_test = df_plp[df_plp['Gene'] == test_gene]['spid'].unique()
        y_true_test = df_valid_test['spid'].isin(carriers_test).astype(int)
        
        result_data['pred_df'] = pd.DataFrame({
            'spid': df_valid_test['spid'],
            f'{test_gene}_score': y_pred_test,
            f'{test_gene}_carrier': y_true_test
        })

    return result_data


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=15, help='Number of parallel processes to use')
    args = parser.parse_args()

    # 1. Load data
    print("Loading training data...")
    df_pheno, df_gene_meta, df_assoc, df_plp, preprocessor = load_training_data(RAW_DATA_PATH, PROCESSED_DATA_DIR)

    df_prev = pd.DataFrame(list(preprocessor.stats_prevalence.items()), columns=['phenotype', 'prevalence'])
    df_mean = pd.DataFrame(list(preprocessor.stats_mean.items()), columns=['phenotype', 'mean'])

    # genes with nominal significant traits > 0
    target_genes = df_gene_meta[df_gene_meta['p_phenotypes'] > 0].index.tolist()
    print(f"Total genes to process: {len(target_genes)}")
    
    # 2. Hyperparameter grid
    pen_candidates = np.round(np.arange(0, 1.05, 0.05), 2)
    conf_candidates = np.round(np.arange(0, 1.05, 0.05), 2)
    param_grid = list(product(pen_candidates, conf_candidates))
    
    # 3. LOG-CV
    max_workers = min(args.workers, os.cpu_count()-1)
    max_workers = os.cpu_count() if max_workers < 1 else max_workers
    print(f"\nStarting parallel LOG-CV with {max_workers} workers for {len(target_genes)} genes...")
    
    results = []
    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=init_worker,
        initargs=(df_pheno, df_assoc, df_plp, df_prev, df_mean)
    ) as executor:
        futures = [executor.submit(process_one_gene_log_cv, g, target_genes, param_grid) for g in target_genes]
        for i, future in enumerate(as_completed(futures)):
            try:
                res = future.result()
                results.append(res)
            except Exception as e:
                print(f"Task failed: {e}")

            if (i+1) % 10 == 0: print(f"Completed {i + 1}/{len(target_genes)} genes...")

    print("Aggregating results...")

    # 1. Best parameters
    param_log = pd.DataFrame([{k:v for k,v in r.items() if k != 'pred_df'} for r in results])
    param_log.to_csv(os.path.join(PROCESSED_DATA_DIR, 'logcv_best_params.csv'), index=False)
    
    # 2. Save predictions
    all_preds = [r['pred_df'].set_index('spid') for r in results if r['pred_df'] is not None]
    if all_preds:
        df_final = pd.concat(all_preds, axis=1)
        df_final = df_final.sort_index()
        output_path = os.path.join(PROCESSED_DATA_DIR, 'gps_predictions_logcv.tsv')
        df_final.to_csv(output_path, sep='\t')
        print(f"Predictions saved to {output_path}")


if __name__ == "__main__":
    main()
