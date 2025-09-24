setwd("/proj/milovelab/wu/scLDAseq")
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)
seurat.file <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/h5ad/sims_", set_level, ".h5seurat")
if (!file.exists(seurat.file)) {
  sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", file_name))
  
  sims <- scuttle::logNormCounts(sims)
  dat <- as.Seurat(sims)
  
  SaveH5Seurat(dat, filename = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/h5ad/sims_", set_level, ".h5seurat"))
}

Convert(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/h5ad/sims_", set_level, ".h5seurat"), dest = "h5ad")
