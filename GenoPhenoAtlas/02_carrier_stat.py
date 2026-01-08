# 02_carrier_stat.py

import os
import pandas as pd


BASE_DIR = "/home/mgbi/projects/gpasd"

# Input file paths
PHENOTYPES_PATH = os.path.join(BASE_DIR, "data/phenotype/spark_clinical_infor_44962.phenotype.csv")
COVARIATES_PATH = os.path.join(BASE_DIR, "data/phenotype/spark_clinical_infor_44962.covariates.csv")
GENOTYPES_PATH = os.path.join(BASE_DIR, "data/genes/output/gene_burden.plp_3726.tsv")

# Output file path
OUTPUT_PATH = os.path.join(BASE_DIR, "data/phewas/output/plp_3726.mut_counts.out.tsv")


def main():
    print("Loading Data...")
    df_phe = pd.read_csv(PHENOTYPES_PATH)
    df_cov = pd.read_csv(COVARIATES_PATH)
    df_gen = pd.read_csv(GENOTYPES_PATH, sep='\t')
    
    id_col = df_phe.columns[0]
    print(f"Detected ID column: {id_col}")

    print("Merging Data...")
    df_merged = df_phe.merge(df_cov, on=id_col).merge(df_gen, on=id_col)
    print(f"Merged DataFrame shape: {df_merged.shape}")

    phe_cols = [c for c in df_phe.columns if c != id_col]
    cov_cols = [c for c in df_cov.columns if c != id_col]
    gen_cols = [c for c in df_gen.columns if c != id_col]

    n_before = len(df_merged)
    df_valid_cov = df_merged.dropna(subset=cov_cols)
    n_after = len(df_valid_cov)
    
    if n_before > n_after:
        print(f"Dropped {n_before - n_after} samples due to missing covariates.")
    
    print("Calculating Gene Mutation Counts...")
    
    results = []

    total_phe = len(phe_cols)
    
    for i, phe in enumerate(phe_cols):
        mask_phe_valid = df_valid_cov[phe].notna()
        df_current_slice = df_valid_cov.loc[mask_phe_valid, gen_cols]

        counts_series = df_current_slice.sum(axis=0)
  
        temp_df = counts_series.reset_index()
        temp_df.columns = ['gene', 'n_mut']
        temp_df['phenotype'] = phe
        
        results.append(temp_df)
        
        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{total_phe} phenotypes...")

    print("Concatenating results...")
    final_result = pd.concat(results, ignore_index=True)

    final_result = final_result[['phenotype', 'gene', 'n_mut']]

    print(f"Saving to {OUTPUT_PATH}...")
    os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_PATH)), exist_ok=True)
    final_result.to_csv(OUTPUT_PATH, sep='\t', index=False)
    
    print("DONE!")


if __name__ == "__main__":
    main()
