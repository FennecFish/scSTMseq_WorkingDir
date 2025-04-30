setwd("/proj/milovelab/wu/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(scater)
library(scran)
library(SC3)
library(ROCR)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims/", file_name))
rowData(sims)$feature_symbol <- rownames(sims)
# remove features with duplicated names
sims <- sims[!duplicated(rowData(sims)$feature_symbol), ]

sims <- quickPerCellQC(sims)
sims <- scuttle::logNormCounts(sims)

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

sims <- sc3(sims, ks = ngroup:ngroup+1, biology = TRUE, n_cores = 3)

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sc3/", 
                                   "sc3_", set_level, ".rds"))
