# src/metrics.py

from typing import List, Dict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from sklearn.metrics import roc_curve, auc, precision_recall_curve

# Random seed for reproducibility
np.random.seed(42)

def calculate_or_and_fisher(
    df: pd.DataFrame, 
    top_percents: List[float], 
    score_col: str = 'score', 
    carrier_col: str = 'carrier',
    correct: bool = False
) -> pd.DataFrame:
    """
    Calculate OR and Fisher's exact test p-value at specified top percentiles.
    """
    results = []
    df_sorted = df.sort_values(score_col, ascending=False)
    total_counts = len(df)
    total_carriers = df[carrier_col].sum()

    for top_percent in top_percents:
        n_top = int(total_counts * top_percent)
        if n_top == 0: continue

        df_top = df_sorted.iloc[:n_top]
        top_counts = len(df_top)
        top_carriers = df_top[carrier_col].sum()
        
        tp = top_carriers
        fp = top_counts - top_carriers
        fn = total_carriers - top_carriers
        tn = (total_counts - total_carriers) - (top_counts - top_carriers)

        table = [[tp, fp], [fn, tn]]

        oddsratio_fisher, p_value = fisher_exact(table, alternative='two-sided')

        if correct:
            tbl = Table2x2(table, shift_zeros=True)
            oddsratio = tbl.oddsratio
            ci_low, ci_upp = tbl.oddsratio_confint()
            is_corrected = 0 in [tp, fp, fn, tn]
        else:
            tbl = Table2x2(table, shift_zeros=False)
            oddsratio = tbl.oddsratio
            try:
                ci_low, ci_upp = tbl.oddsratio_confint()
            except:
                ci_low, ci_upp = np.nan, np.nan
            is_corrected = False

        top_yield = top_carriers / top_counts if top_counts > 0 else np.nan
        total_yield = total_carriers / total_counts if total_counts > 0 else np.nan

        results.append({
            'top_percent': f"{top_percent:.1%}",
            'n_top': top_counts,
            'n_top_carriers': top_carriers,
            'top_yield': top_yield,
            'n_total': total_counts,
            'n_total_carriers': total_carriers,
            'total_yield': total_yield,
            'or': oddsratio,
            'or_ci_low': ci_low,
            'or_ci_upp': ci_upp,
            'p': p_value,
            'is_corrected': is_corrected
        })

    return pd.DataFrame(results)


def bootstrap_auc(
    y_true: np.ndarray, 
    y_scores: np.ndarray, 
    n_bootstrap: int = 1000,
    apr: bool = False
) -> Dict:

    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)
    n_samples = len(y_true)
    
    auc_scores, prauc_scores = [], []
    tprs, precisions = [], []
    base_fpr = np.linspace(0, 1, 101)
    base_recall = np.linspace(0, 1, 101)
    
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    
    if len(pos_idx) == 0 or len(neg_idx) == 0:
        return {}

    for _ in range(n_bootstrap):
        # Bootstrap
        indices = np.random.choice(n_samples, n_samples, replace=True)

        y_true_bs = y_true[indices]
        y_score_bs = y_scores[indices]
        
        if len(np.unique(y_true_bs)) < 2: continue
        
        if apr:
            precision, recall, _ = precision_recall_curve(y_true_bs, y_score_bs)
            prauc = auc(recall, precision)
            prauc_scores.append(prauc)
            prec_interp = np.interp(base_recall, recall[::-1], precision[::-1])
            precisions.append(prec_interp)

        fpr, tpr, _ = roc_curve(y_true_bs, y_score_bs)
        auc_score = auc(fpr, tpr)
        auc_scores.append(auc_score)
        tpr_interp = np.interp(base_fpr, fpr, tpr)
        tprs.append(tpr_interp)

    tprs = np.array(tprs)
    precisions = np.array(precisions)

    res_dict = {
        'auc_mean': np.mean(auc_scores),
        'auc_median': np.median(auc_scores),
        'auc_ci_lower': np.percentile(auc_scores, 2.5),
        'auc_ci_upper': np.percentile(auc_scores, 97.5),
        'base_fpr': base_fpr.tolist(),
        'tpr_mean': np.mean(tprs, axis=0).tolist(),
        'tpr_lower': np.percentile(tprs, 2.5, axis=0).tolist(),
        'tpr_upper': np.percentile(tprs, 97.5, axis=0).tolist()
    }
    if apr:
        res_dict.update({
            'prauc_mean': np.mean(prauc_scores),
            'prauc_median': np.median(prauc_scores),
            'prauc_ci_lower': np.percentile(prauc_scores, 2.5),
            'prauc_ci_upper': np.percentile(prauc_scores, 97.5),
            'base_recall': base_recall.tolist(),
            'precision_mean': np.mean(precisions, axis=0).tolist(),
            'precision_lower': np.percentile(precisions, 2.5, axis=0).tolist(),
            'precision_upper': np.percentile(precisions, 97.5, axis=0).tolist()
        })

    return res_dict
