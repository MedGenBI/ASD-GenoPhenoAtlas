# scripts/03_gene_level_eval.py

import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.metrics import bootstrap_auc, calculate_or_and_fisher
from src.config import PROCESSED_DATA_DIR

def main():
    # 1. Load predictions
    pred_file = os.path.join(PROCESSED_DATA_DIR, 'gps_predictions_logcv.tsv')
    if not os.path.exists(pred_file):
        raise FileNotFoundError("Prediction file not found. Run 02_model_optimization.py first.")
    
    print("Loading predictions...")
    df_preds = pd.read_csv(pred_file, sep='\t', index_col='spid')
    
    # 2. Extract genes: {gene}_score, {gene}_carrier
    genes = set([col.replace('_score', '') for col in df_preds.columns if col.endswith('_score')])
    print(f"\nEvaluating {len(genes)} genes...")
    
    gene_stats = []
    
    for i, gene in enumerate(genes):
        score_col = f"{gene}_score"
        carrier_col = f"{gene}_carrier"

        subset = df_preds[[score_col, carrier_col]].dropna()
        if len(subset) == 0: continue
        
        y_true = subset[carrier_col].values
        y_score = subset[score_col].values
        
        if y_true.sum() == 0: continue
        
        # Calculate AUC, APR
        auc_res = bootstrap_auc(y_true, y_score, n_bootstrap=1000, apr=True)

        scalar_res = {k: v for k, v in auc_res.items() if not isinstance(v, (np.ndarray, list))}
        
        # Calculate top 5% odds ratio
        top_percents = [0.05]
        df_top_or = calculate_or_and_fisher(subset, top_percents, score_col=score_col, carrier_col=carrier_col)
        or_dict = df_top_or.iloc[0].to_dict()

        record = {
            'gene': gene,
            'n_samples': len(subset),
            'n_carriers': y_true.sum(),
            **scalar_res,
            **or_dict
        }
        gene_stats.append(record)
        
        if (i+1) % 10 == 0:
            print(f"Evaluated {i+1}/{len(genes)} genes...")

    df_stats = pd.DataFrame(gene_stats)

    df_gene_meta = pd.read_csv(os.path.join(PROCESSED_DATA_DIR, 'gene_metadata.csv'))
    df_stats = df_stats.merge(df_gene_meta[['gene', 'p_phenotypes', 'fdr_phenotypes']], on='gene', how='left')
    df_stats = df_stats.sort_values(by=['fdr_phenotypes', 'p_phenotypes', 'gene'], ascending=[False, False, True])

    # 3. Save results
    output_path = os.path.join(PROCESSED_DATA_DIR, 'gene_level_evaluation_metrics.tsv')
    df_stats.to_csv(output_path, sep='\t', index=False)
    print(f"Gene-level metrics saved to {output_path}")


if __name__ == "__main__":
    main()
