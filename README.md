# ASD-GenoPhenoAtlas

## Overview

This repository contains the source code and analysis scripts for the study **"Genotype-phenotype atlas dissects autism heterogeneity and informs gene prioritization"**. It provides a comprehensive framework for constructing a genotype-phenotype association atlas in Autism Spectrum Disorder (ASD) and implements the Gene-centric Polyphenotypic Score (GPS), a novel quantitative metric that aggregates phenotype signals to prioritize genes for P/LP variant carriers.

## Key Components

- **Association analysis**: Gene-based genotype-phenotype association analysis and post-processing.
- **GPS**: Full implementation of the GPS, including parameter optimization and performance evaluation.
- **Visualization**: Scripts used to generate or support the main figures presented in the manuscript.

## Installation

Clone this repo:

```bash
git clone https://github.com/MedGenBI/ASD-GenoPhenoAtlas.git
cd ASD-GenoPhenoAtlas
```

Create a Conda environment:

```bash
conda env create -f environment.yml
conda activate gpasd
```

## Reproducing the analyses

If you have authorized access to the required SPARK data, prepare the study inputs as described below and follow the corresponding workflow steps.

For reproducibility, the random seed is set to 42 by default across all relevant workflow steps.

### Data preparation

Input files that can be shared are provided directly in `Data/`. For participant-level inputs subject to data-use restrictions, format examples are provided in `Data/examples/`; authorized users should prepare the corresponding files following these examples, retain the same filenames, and place them in `Data/` before running the analyses.

The example files in `Data/examples/` have the following formats:

- `phenotype.csv`: SPARK ES phenotype matrix containing a subject ID column (`spid`) and multiple trait columns. The trait column names correspond to the `trait_id`s defined in `trait_metadata.csv`. Binary traits are coded as 0, 1, or missing, while quantitative traits retain their original numeric values.

- `phenotype_GS.csv`: SPARK GS phenotype matrix, sharing the identical column structure (`spid` and trait columns) as the ES matrix.

- `phenotype_binary_coded.csv`: SPARK ES binary severity matrix, encoded with 0, 1, or missing values. The conversion mapping from original values to binary indicators is detailed in the manuscript Methods.

- `phenotype_binary_coded_GS.csv`: SPARK GS binary severity matrix, derived using the exact same coding scheme as `phenotype_binary_coded.csv`.

- `covariate.csv`: Covariate matrix comprising the subject ID column (`spid`) and the corresponding covariate variables used for adjustment.

For the GPS component, input paths within the shared `Data/` directory can be configured in two ways: edit the defaults in `GPS/src/config.py`, or set environment variables. Environment variables override the corresponding values in `config.py`. To configure these paths using environment variables, you can use the following commands:

```bash
export GPS_DATA_DIR=/path/to/ASD-GenoPhenoAtlas/Data
export GPS_TRAIT_METADATA=/path/to/ASD-GenoPhenoAtlas/Data/trait_metadata.csv
export GPS_PHENOTYPE_PATH=/path/to/ASD-GenoPhenoAtlas/Data/phenotype.csv
export GPS_PLP_PATH=/path/to/ASD-GenoPhenoAtlas/Data/ACMG_PLP_variants.tsv
export GPS_WGS_PHENO_PATH=/path/to/ASD-GenoPhenoAtlas/Data/phenotype_GS.csv
export GPS_WGS_PLP_PATH=/path/to/ASD-GenoPhenoAtlas/Data/ACMG_PLP_variants_GS.tsv
```

### Reproducing the genotype-phenotype association analysis

Enter the directory:

```bash
cd /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas
```

1. Run association analysis:

```bash
python gpasd.py association \
  --pheno /path/to/ASD-GenoPhenoAtlas/Data/phenotype.csv \
  --geno /path/to/ASD-GenoPhenoAtlas/Data/burden_genotype.tsv \
  --cov /path/to/ASD-GenoPhenoAtlas/Data/covariate.csv \
  --trait-metadata /path/to/ASD-GenoPhenoAtlas/Data/trait_metadata.csv \
  --cores 16 \
  --output /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.raw.tsv
```

`--pheno` specifies the phenotype matrix, `--geno` specifies the burden-genotype matrix, `--cov` specifies the covariate matrix, and `--trait-metadata` specifies the phenotype metadata file. `--cores` sets the number of parallel processes, and `--output` specifies the raw association-results TSV file.

2. Generate permuted genotype files:

```bash
python gpasd.py geno-permute \
  --geno /path/to/ASD-GenoPhenoAtlas/Data/burden_genotype.tsv \
  --outdir /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/permutation/geno \
  --n-perms 1000 \
  --seed 42
```

`--geno` specifies the original burden-genotype file, and `--outdir` specifies the directory for the generated null genotype files. The default is 1,000 permutations with random seed 42.

3. Run association analysis for the permuted genotypes:

```bash
python gpasd.py association \
  --pheno /path/to/ASD-GenoPhenoAtlas/Data/phenotype.csv \
  --cov /path/to/ASD-GenoPhenoAtlas/Data/covariate.csv \
  --trait-metadata /path/to/ASD-GenoPhenoAtlas/Data/trait_metadata.csv \
  --geno-dir /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/permutation/geno \
  --outdir /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/permutation/assoc \
  --cores 48
```

`--pheno` specifies the phenotype matrix, `--cov` specifies the covariate matrix, and `--trait-metadata` specifies the phenotype metadata file. `--geno-dir` specifies the directory containing the permuted genotype files, `--outdir` specifies the output directory, and `--cores` sets the number of parallel processes.

4. Apply multiple-testing correction:

```bash
python gpasd.py multitest \
  --assoc-file /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.raw.tsv \
  --null-pvals-dir /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/permutation/assoc \
  --trait-metadata /path/to/ASD-GenoPhenoAtlas/Data/trait_metadata.csv \
  --output /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.multitest.tsv
```

`--assoc-file` specifies the raw association-results file, `--null-pvals-dir` specifies the directory containing the permutation association results, `--trait-metadata` specifies the phenotype metadata file, and `--output` specifies the multiple-testing results file.

5. Format the results:

```bash
python gpasd.py format \
  -i /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.multitest.tsv \
  -o /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.tsv \
  --trait-metadata /path/to/ASD-GenoPhenoAtlas/Data/trait_metadata.csv \
  --sort
```

`-i` specifies the multiple-testing results file, `-o` specifies the formatted output file, `--trait-metadata` specifies the phenotype metadata file, and `--sort` sorts the formatted results by P value.


### Reproducing the GPS workflow

Enter the GPS directory:

```bash
cd /path/to/ASD-GenoPhenoAtlas/GPS
```

1. Prepare the derivation data:

```bash
python scripts/01_data_preparation.py \
  --association-file /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.tsv \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

`--association-file` specifies the genotype-phenotype association results file, and `--outdir` specifies the output directory of GPS.

2. Optimize the GPS model and generate derivation predictions:

```bash
python scripts/02_model_optimization.py \
  --association-file /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.tsv \
  --workers 16 \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

`--association-file` specifies the genotype-phenotype association results file, `--workers` sets the number of parallel workers, and `--outdir` specifies the output directory of GPS.

3. Evaluate gene-level discrimination:

```bash
python scripts/03_gene_level_eval.py \
  --workers 16 \
  --n-bootstrap 1000 \
  --seed 42 \
  --pred-file /path/to/ASD-GenoPhenoAtlas/GPS/output/02_model_optimization/output/gps_predictions_logcv.tsv \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

`--pred-file` specifies the step-02 prediction table, `--workers` sets the number of parallel workers, `--n-bootstrap` sets the bootstrap resamples per gene, `--seed` sets the base random seed, and `--outdir` specifies the output directory of GPS.

4. Compare GPS with the severity composite (phenotype burden score):

```bash
python scripts/04_gps_vs_severity.py \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output \
  --phenotype-file /path/to/ASD-GenoPhenoAtlas/Data/phenotype_binary_coded.csv \
  --pred-file /path/to/ASD-GenoPhenoAtlas/GPS/output/02_model_optimization/output/gps_predictions_logcv.tsv \
  --gene-level-metrics-file /path/to/ASD-GenoPhenoAtlas/GPS/output/03_gene_level_eval/output/gene_level_evaluation_metrics.tsv \
  --formats png,pdf,svg
```

`--phenotype-file` specifies the binary severity matrix, `--pred-file` specifies the step-02 GPS predictions, `--gene-level-metrics-file` specifies the step-03 metrics, `--formats` selects the figure formats, and `--outdir` specifies the output directory of GPS.

5. Calculate gene-level carrier enrichment and create the candidate panel:

```bash
python scripts/05_gene_level_or.py \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output \
  --pred-file /path/to/ASD-GenoPhenoAtlas/GPS/output/02_model_optimization/output/gps_predictions_logcv.tsv \
  --gene-metrics-file /path/to/ASD-GenoPhenoAtlas/GPS/output/03_gene_level_eval/output/gene_level_evaluation_metrics.tsv \
  --candidate-genes auto \
  --top-percents 0.2,0.1,0.05,0.01
```

`--pred-file` specifies the step-02 prediction table, `--gene-metrics-file` specifies the step-03 metrics, `--top-percents` sets the top fractions, and `--outdir` specifies the output directory of GPS.

6. Generate derivation-cohort patient-level predictions:

```bash
python scripts/06_patient_level_pred_drv.py \
  --candidate-genes /path/to/ASD-GenoPhenoAtlas/GPS/output/05_gene_level_or/output/selected_genes.txt \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

The command uses the candidate gene panel selected in step-05 and the shared output directory.

7. Generate external-cohort patient-level predictions:

```bash
python scripts/07_patient_level_pred_ext.py \
  --association-file /path/to/ASD-GenoPhenoAtlas/GenoPhenoAtlas/output/geno_pheno_association_results.tsv \
  --candidate-genes /path/to/ASD-GenoPhenoAtlas/GPS/output/05_gene_level_or/output/selected_genes.txt \
  --workers 16 \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

`--association-file` specifies the genotype-phenotype association results file, `--candidate-genes` specifies the panel from step-05, and `--workers` sets the number of parallel workers.

8. Evaluate derivation- and external-cohort patient-level predictions against benchmark results:

```bash
python scripts/08_patient_level_eval.py \
  --benchmark-summary /path/to/ASD-GenoPhenoAtlas/Data/benchmark_drv.tsv \
  --gps-prediction /path/to/ASD-GenoPhenoAtlas/GPS/output/06_patient_level_pred_drv/output/patient_level_predication.tsv \
  --analysis-label derivation \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

```bash
python scripts/08_patient_level_eval.py \
  --benchmark-summary /path/to/ASD-GenoPhenoAtlas/Data/benchmark_ext.tsv \
  --gps-prediction /path/to/ASD-GenoPhenoAtlas/GPS/output/07_patient_level_pred_ext/output/patient_level_predication.tsv \
  --analysis-label external \
  --outdir /path/to/ASD-GenoPhenoAtlas/GPS/output
```

`--benchmark-summary` specifies the non-GPS methods benchmark summary, `--gps-prediction` specifies the corresponding GPS prediction table, `--analysis-label` identifies the derivation or external cohort, and `--outdir` specifies the output directory of GPS.

## Data

All raw genetic and phenotypic data for the SPARK cohorts are available upon request from [SFARI Base](https://base.sfari.org/). The full variant discovery, quality control, annotation, polygenic score (PGS) calculation, and other related procedures were conducted using publicly available software and datasets as fully detailed in the manuscript.

## Citation

Citation information will be updated upon publication.
