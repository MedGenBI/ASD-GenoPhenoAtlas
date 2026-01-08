# src/dataloader.py

import pandas as pd
import pickle
import os
from typing import Tuple
from .preprocessor import PhenotypePreprocessor

def load_raw_data(association_path: str, phenotype_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:

    print(f"Loading association data from: {association_path}")
    assoc_df = pd.read_excel(association_path, sheet_name='PLP Association')
    
    print(f"Loading phenotype data from: {phenotype_path}")
    pheno_df = pd.read_csv(phenotype_path)
    
    return assoc_df, pheno_df

def generate_metadata(assoc_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

    print("Generating metadata statistics (gene_count, phe_count)...")

    # 1. Gene counts
    gene_stats = assoc_df.groupby('gene').agg(
        p_phenotypes=('p', lambda x: (x < 0.05).sum()),
        fdr_phenotypes=('p_fdr', lambda x: (x < 0.05).sum())
    )
    gene_count_df = gene_stats.sort_values(by=['fdr_phenotypes', 'p_phenotypes'], ascending=False)

    # 2. Phenotype counts
    phe_stats = assoc_df.groupby('phenotype').agg(
        p_genes=('p', lambda x: (x < 0.05).sum()),
        fdr_genes=('p_fdr', lambda x: (x < 0.05).sum())
    )
    phe_count_df = phe_stats.sort_values(by='fdr_genes', ascending=False)

    return gene_count_df, phe_count_df

def load_training_data(
    raw_data_path: str, 
    processed_data_dir: str
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, PhenotypePreprocessor]:

    # 1. Load processed phenotype data
    pheno_path = os.path.join(processed_data_dir, 'phenotypes_processed_int.csv')
    df_pheno = pd.read_csv(pheno_path)
    
    # 2. Metadata
    gene_meta_path = os.path.join(processed_data_dir, 'gene_metadata.csv')
    df_gene_meta = pd.read_csv(gene_meta_path, index_col='gene')
    
    # 3. Association data
    df_assoc = pd.read_excel(raw_data_path, sheet_name='PLP Association')
    
    # 4. PLP carrier info
    df_plp = pd.read_excel(raw_data_path, sheet_name='PLP ACMG')
    
    # 5. Preprocessor
    model_path = os.path.join(processed_data_dir, 'preprocessor_model.pkl')
    with open(model_path, 'rb') as f:
        preprocessor = pickle.load(f)
        
    print("All training data loaded successfully.")
    return df_pheno, df_gene_meta, df_assoc, df_plp, preprocessor
