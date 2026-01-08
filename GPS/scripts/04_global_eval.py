# scripts/04_global_eval.py

import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.metrics import bootstrap_auc, calculate_or_and_fisher
from src.plotting import plot_global_roc
from src.config import PROCESSED_DATA_DIR, OUTPUT_DIR

np.random.seed(42)

def create_global_rank_df(df_preds, gene_list):

    dfs = []
    for gene in gene_list:
        score_col = f"{gene}_score"
        carrier_col = f"{gene}_carrier"
        if score_col not in df_preds.columns: continue
        
        sub_df = df_preds[[score_col, carrier_col]].dropna().copy()
        sub_df.columns = ['score', 'carrier']
        sub_df['gene'] = gene
        sub_df['rank'] = sub_df['score'].rank(pct=True)
        sub_df['spid'] = sub_df.index
        
        dfs.append(sub_df)
    
    if not dfs:
        return pd.DataFrame()
    
    return pd.concat(dfs, ignore_index=True)


def main():

    print("Global Evaluation...")

    # 1. Load data
    pred_file = os.path.join(PROCESSED_DATA_DIR, 'gps_predictions_logcv.tsv')
    metric_file = os.path.join(PROCESSED_DATA_DIR, 'gene_level_evaluation_metrics.tsv')
    
    df_preds = pd.read_csv(pred_file, sep='\t', index_col='spid')
    df_metrics = pd.read_csv(metric_file, sep='\t')
    

    # 2. Grouping
    # High-confidence: AUC > 0.8 & FDR Phenotypes > 0
    hc_genes = df_metrics[
        (df_metrics['auc_median'] >= 0.8) & 
        (df_metrics['fdr_phenotypes'] > 0)
    ]['gene'].tolist()
    
    # Extended: AUC > 0.8
    ex_genes = df_metrics[df_metrics['auc_median'] >= 0.8]['gene'].tolist()
    
    print(f"High-confidence Set: {len(hc_genes)} genes")
    print(f"Extended Set: {len(ex_genes)} genes")
    
    with open(os.path.join(PROCESSED_DATA_DIR, 'hc_gene_list.txt'), 'w') as f_hc:
        for gene in hc_genes:
            f_hc.write(f"{gene}\n")
    with open(os.path.join(PROCESSED_DATA_DIR, 'ex_gene_list.txt'), 'w') as f_ex:
        for gene in ex_genes:
            f_ex.write(f"{gene}\n")


    # 3. Create global rank data
    df_rank_hc = create_global_rank_df(df_preds, hc_genes)
    df_rank_ex = create_global_rank_df(df_preds, ex_genes)
    

    # 4. Save rank data
    df_rank_hc.to_csv(os.path.join(PROCESSED_DATA_DIR, 'wes_global_rank_hc.tsv'), sep='\t', index=False)
    df_rank_ex.to_csv(os.path.join(PROCESSED_DATA_DIR, 'wes_global_rank_ex.tsv'), sep='\t', index=False)
    

    # 5. Global Evaluation (AUC and APR)
    print("\nCalculating Global AUC & APR...")

    global_auc_hc = bootstrap_auc(df_rank_hc['carrier'], df_rank_hc['rank'], apr=True)
    global_auc_ex = bootstrap_auc(df_rank_ex['carrier'], df_rank_ex['rank'], apr=True)
    
    # Save results
    pd.DataFrame([global_auc_hc]).to_csv(os.path.join(PROCESSED_DATA_DIR, 'wes_global_auc_hc.csv'), index=False)
    pd.DataFrame([global_auc_ex]).to_csv(os.path.join(PROCESSED_DATA_DIR, 'wes_global_auc_ex.csv'), index=False)
    
    # Plotting
    plot_path = os.path.join(OUTPUT_DIR, 'plots', 'fig4.c.AUC.wes.pdf')
    plot_global_roc(global_auc_hc, global_auc_ex, plot_path)
    

    # 6. Global Evaluation (OR)
    print("\nCalculating Odds Ratios...")
    top_percents = [0.2, 0.1, 0.05, 0.01, 0.001]
    or_hc = calculate_or_and_fisher(df_rank_hc, top_percents, score_col='rank', carrier_col='carrier')
    or_hc['group'] = 'High-confidence'
    
    or_ex = calculate_or_and_fisher(df_rank_ex, top_percents, score_col='rank', carrier_col='carrier')
    or_ex['group'] = 'Extended'
    
    df_or_all = pd.concat([or_hc, or_ex])
    or_out_path = os.path.join(OUTPUT_DIR, 'plots', 'data', 'fig4.c.wes.global_or.forest.tsv')
    os.makedirs(os.path.dirname(or_out_path), exist_ok=True)
    df_or_all.to_csv(or_out_path, index=False)
    print(f"Global ORs saved to {or_out_path}")


if __name__ == "__main__":
    main()
