# scripts/01_data_preparation.py

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.append(project_root)

from src.config import (
    RAW_DATA_PATH, PHENOTYPE_PATH, PROCESSED_DATA_DIR,
    CONTINUOUS_TRAITS, BINARY_TRAITS
)
from src.dataloader import load_raw_data, generate_metadata
from src.preprocessor import PhenotypePreprocessor

def main():
    # 1. Data loading
    print("Loading Data...")
    assoc_df, pheno_df = load_raw_data(RAW_DATA_PATH, PHENOTYPE_PATH)
    
    # 2. Basic data construction
    print("\nGenerating Metadata (Gene & Phenotype Counts)...")
    gene_count_df, phe_count_df = generate_metadata(assoc_df)
    
    gene_count_df.to_csv(os.path.join(PROCESSED_DATA_DIR, 'gene_metadata.csv'))
    phe_count_df.to_csv(os.path.join(PROCESSED_DATA_DIR, 'phenotype_metadata.csv'))
    print("Metadata saved.")

    # 3. Phenotype preprocessing
    print("\nPreprocessing Phenotypes...")
    
    preprocessor = PhenotypePreprocessor(
        continuous_traits=CONTINUOUS_TRAITS,
        binary_traits=BINARY_TRAITS,
        missing_threshold=0.2
    )
    pheno_processed = preprocessor.fit_transform(pheno_df)
    
    print(f"Processed phenotype data shape: {pheno_processed.shape}")
    
    # 4. Save
    output_pheno_path = os.path.join(PROCESSED_DATA_DIR, 'phenotypes_processed_int.csv')
    pheno_processed.to_csv(output_pheno_path, index=False)
    print(f"Processed phenotype data saved to: {output_pheno_path}")

    model_path = os.path.join(PROCESSED_DATA_DIR, 'preprocessor_model.pkl')
    preprocessor.save(model_path)
    
    # Preview
    print("\n--- Statistics Preview ---")
    print(f"Example Continuous Mean (INT): {list(preprocessor.stats_mean.items())[:3]}")
    print(f"Example Binary Prevalence: {list(preprocessor.stats_prevalence.items())[:3]}")

if __name__ == "__main__":
    main()
