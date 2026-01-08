# scripts/06_external_validation_eval.py

import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.config import OUTPUT_DIR
from src.metrics import bootstrap_auc, calculate_or_and_fisher
from src.plotting import plot_global_roc

np.random.seed(42)

def process_rank_conversion(df_pred, label="Set"):
    dfs = []
    genes = [c.replace('_score', '') for c in df_pred.columns if c.endswith('_score')]
    
    sample_count = 0
    gene_count_valid = 0
    
    for gene in genes:
        score_col = f"{gene}_score"
        carrier_col = f"{gene}_carrier"
        
        subset = df_pred[[score_col, carrier_col]].dropna().copy()
        if subset.empty: continue
        
        subset.columns = ['score', 'carrier']
        subset['rank'] = subset['score'].rank(pct=True)
        subset['gene'] = gene
        subset['spid'] = subset.index
        dfs.append(subset)
        
        sample_count += subset['carrier'].sum()
        gene_count_valid += 1
        
    print(f"({label}) Final PLP Carriers: {int(sample_count)}; Genes Involved: {gene_count_valid}")
    
    if not dfs: return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def main():
    print("External Validation - Evaluation...")
    
    # 1. Load predictions
    hc_path = os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_predictions_hc.tsv')
    ex_path = os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_predictions_ex.tsv')
    
    df_hc = pd.read_csv(hc_path, sep='\t', index_col='spid')
    df_ex = pd.read_csv(ex_path, sep='\t', index_col='spid')
    
    # 2. Create global rank 
    df_rank_hc = process_rank_conversion(df_hc, "High-confidence")
    df_rank_ex = process_rank_conversion(df_ex, "Extended")

    df_rank_hc.to_csv(os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_global_rank_hc.tsv'), sep='\t', index=False)
    df_rank_ex.to_csv(os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_global_rank_ex.tsv'), sep='\t', index=False)
    
    # 3. Calculate global AUC & APR
    print("\nCalculating AUC & APR...")
    res_hc = bootstrap_auc(df_rank_hc['carrier'], df_rank_hc['rank'], stratify=False, pr=True)
    res_ex = bootstrap_auc(df_rank_ex['carrier'], df_rank_ex['rank'], stratify=False, pr=True)
    
    # Save results
    pd.DataFrame([res_hc]).to_csv(os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_global_auc_hc.csv'))
    pd.DataFrame([res_ex]).to_csv(os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_global_auc_ex.csv'))
    
    # 4. Plotting
    fig_path = os.path.join(OUTPUT_DIR, 'plots', 'fig4.d.AUC.wgs.pdf')
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plot_global_roc(res_hc, res_ex, fig_path)
    
    # 5. Calculate OR
    print("\nCalculating OR...")
    top_percents = [0.2, 0.1, 0.05, 0.01, 0.001]
    
    or_hc = calculate_or_and_fisher(df_rank_hc, top_percents, score_col='rank', carrier_col='carrier')
    or_hc['group'] = 'High-confidence'
    
    or_ex = calculate_or_and_fisher(df_rank_ex, top_percents, score_col='rank', carrier_col='carrier')
    or_ex['group'] = 'Extended'
    
    df_or_all = pd.concat([or_hc, or_ex])
    df_or_all.to_csv(os.path.join(OUTPUT_DIR, 'processed_data', 'wgs_global_or.tsv'), sep='\t', index=False)
    print("WGS Evaluation completed.")


if __name__ == "__main__":
    main()
