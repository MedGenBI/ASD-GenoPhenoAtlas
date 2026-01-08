# 01_geno_pheno_association_analysis.R

#### R package
library(PheWAS)

base_dir <- "/home/mgbi/projects/gpasd"

# Set working directory
setwd(base_dir)

# Path settings
phenotypes_path <- file.path(base_dir, "data/phenotype/spark_clinical_infor_44962.phenotype.csv")
covariates_path <- file.path(base_dir, "data/phenotype/spark_clinical_infor_44962.covariates.csv")

genotypes_path <- 'data/genes/output/gene_burden.plp_3726.tsv'

additive.genotypes <- F
pheno.transform <- 'int'
cores <- 16

use.logistf <- F

if (use.logistf) {
    phewas_results <- 'data/phewas/output/plp_3726.logistf.phewas.out.tsv'
} else {
    phewas_results <- 'data/phewas/output/plp_3726.phewas.out.tsv'
}


#### Phenotypes
phenotypes <- read.table(phenotypes_path, header = T, sep = ",", check.names = F)

## Binary traits
for (i in 32:79) {
  phenotypes[i] <- as.logical(phenotypes[[i]])
}

## Continuous traits
pheno.continuous <- colnames(phenotypes)[2:31]
if (pheno.transform == 'int') {
  ## rank-based inverse-normal transformation
  message("perform rank-based inverse-normal transformation for continuous variables!")
  for (p.con in pheno.continuous) {
    x <- phenotypes[p.con]
    phenotypes[p.con] <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  }
} else if (pheno.transform == 'log') {
  ## log transformation
  message("perform log transformation for continuous variables!")
  for (p.con in pheno.continuous) {
    x <- phenotypes[p.con]
    phenotypes[p.con] <- log(1+x)
  }
} else if (pheno.transform == 'zscore') {
  ## z-score transformation
  message("perform z-score normalization for continuous variables!")
  for (p.con in pheno.continuous) {
    x <- phenotypes[[p.con]]
    phenotypes[p.con] <- (x-mean(na.omit(x))) / sd(na.omit(x))
  }
}

#### Covariates
covariates <- read.table(covariates_path, header = T, sep = ",", check.names = F)

for (i in c(2,4,5)) {
  covariates[i] <- as.factor(covariates[[i]])
}

#### Genotypes
genotypes <- read.table(genotypes_path, header = T, sep = "\t", check.names = F)

#### Analysis
if (use.logistf) {
  
  id=intersect(names(phenotypes),names(genotypes))
  data = merge(phenotypes,genotypes,by=id)
  data = merge(data,covariates,by=id)
  
  phenotypes_cols = names(phenotypes)
  phenotypes_cols = phenotypes_cols[!(phenotypes_cols %in% id)]
  genotypes_cols = names(genotypes)
  genotypes_cols = genotypes_cols[!(genotypes_cols %in% id)]
  covariates_cols = names(covariates)
  covariates_cols = covariates_cols[!(covariates_cols %in% id)]
  
  if(class(covariates_cols)!="list") { covariates_cols=list(covariates_cols)}
  
  full_list=data.frame(t(expand.grid(phenotypes_cols,genotypes_cols,covariates_cols,stringsAsFactors=F)),stringsAsFactors=F)
  
  results = phewas_ext(
    phenotypes_cols, genotypes_cols, data, covariates = covariates_cols, cores = cores, 
    additive.genotypes = additive.genotypes, method="logistf"
  )
} else {
  results = phewas(
    phenotypes, genotypes, covariates = covariates, cores = cores, 
    additive.genotypes = additive.genotypes
  )
}

#### Save results
write.table(results, phewas_results, sep = "\t", row.names = F)
