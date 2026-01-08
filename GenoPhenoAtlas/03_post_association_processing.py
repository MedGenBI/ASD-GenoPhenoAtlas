# 03_post_association_processing.py

import os
import pandas as pd
from scipy import stats


BASE_DIR = "/home/mgbi/projects/gpasd"

ASSOC_PATH = os.path.join(BASE_DIR, f'data/phewas/output/plp_3726.phewas.out.tsv')
LOGISTF_PATH = os.path.join(BASE_DIR, f'data/phewas/output/plp_3726.logistf.phewas.out.tsv')
MUT_COUNTS_PATH = os.path.join(BASE_DIR, f'data/phewas/output/plp_3726.mut_counts.out.tsv')

OUTPUT_PATH = os.path.join(BASE_DIR, f'data/phewas/output/plp_3726.phewas.processed.tsv')

MUT_COUNT_CUTOFF = 5
FIRTH_P_THRESHOLD = 0.01

COLS_TO_USE = ['phenotype', 'snp', 'beta', 'SE', 'OR', 'p', 'type', 'n_total', 'n_cases', 'n_controls']

def load_and_clean_data(filepath, cols):

    print(f"Loading {filepath}...")
    df = pd.read_csv(filepath, sep='\t', usecols=cols)
    original_len = len(df)
    df = df.dropna(subset=['p']).copy()

    # snp -> gene
    df['snp'] = df['snp'].str.strip("`")
    df.rename(columns={"snp": "gene"}, inplace=True)

    print(f"  - Loaded {len(df)} rows (dropped {original_len - len(df)} rows with NaN p-value)")
    return df

def main():

    print("Loading Data...")
    df_assoc = load_and_clean_data(ASSOC_PATH, COLS_TO_USE)
    df_logistf = load_and_clean_data(LOGISTF_PATH, COLS_TO_USE)
    

    print(f"Loading mutation counts from {MUT_COUNTS_PATH}...")
    df_mut = pd.read_csv(MUT_COUNTS_PATH, sep="\t")
    df_merged = pd.merge(df_assoc, df_mut, on=["phenotype", "gene"], how="left")
        
    print(f"  - Merged shape: {df_merged.shape}")

    n_before = len(df_merged)
    df_merged = df_merged[df_merged['n_mut'] >= MUT_COUNT_CUTOFF].copy()
    print(f"  - Filtered {n_before - len(df_merged)} pairs with n_mut < {MUT_COUNT_CUTOFF}. Final shape: {df_merged.shape}")


    print("Applying Firth Regression Correction (p < 0.01)...")
    df_merged.set_index(['phenotype', 'gene'], inplace=True)
    df_logistf.set_index(['phenotype', 'gene'], inplace=True)

    mask_correction = df_merged['p'] < FIRTH_P_THRESHOLD
    idx_correction = df_merged[mask_correction].index
    idx_valid_correction = idx_correction.intersection(df_logistf.index)
    
    if len(idx_valid_correction) > 0:

        cols_to_update = ['beta', 'SE', 'OR', 'p', 'type']
        df_merged.update(df_logistf.loc[idx_valid_correction, cols_to_update])
        print(f"  - Updated {len(idx_valid_correction)} rows with Firth regression results.")
    else:
        print("  - No rows required Firth correction.")
        
    df_merged.reset_index(inplace=True)


    print("Performing FDR Correction (Benjamini-Hochberg)...")
    df_merged = df_merged.dropna(subset=['p']).copy()
    
    df_merged['p_fdr'] = stats.false_discovery_control(df_merged['p'], method='bh')
    
    final_cols = ['phenotype', 'gene', 'beta', 'SE', 'OR', 'p', 'p_fdr', 'n_mut', 'type', 'n_total', 'n_cases', 'n_controls']
    final_cols = [c for c in final_cols if c in df_merged.columns]
    df_final = df_merged[final_cols]


    print(f"Saving results to {OUTPUT_PATH}...")
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    df_final.to_csv(OUTPUT_PATH, sep='\t', index=False)
    print("DONE!")


if __name__ == "__main__":
    main()
