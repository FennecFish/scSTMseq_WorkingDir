# This script is to run propeller with scSTMseq output
# with replicates as block in limma

setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
library(speckle)
library(limma)
library(CellBench)
library(mclust)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
scSTM_dir <- args[2]
file_name <- args[3]
scSTM_name <- file_name

file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat("scSTM dir is",scSTM_dir, "\n")
cat(file_name, "\n")

set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)

################################################################################
############################## read in seurat ##################################
################################################################################
# select ones with the best clustering results
scSTMobj <- readRDS(scSTM_name)
source("R/useful_functions.R")
source("doc/simulation/1MultiSample_Simulation/benchmark/functions_replicates.R")
# select the scSTMobj with the largest bound
scSTMobj <- select_top_scSTM(scSTMobj)
# sample from the posterior distribution
thetaRep <- PosteriorReplicates(model = scSTMobj, nsims = 100)

propRep.logit <- lapply(thetaRep, function(x){
  TransformReplicates(x, "logit")
})
names(propRep.logit) <- paste0("Replicates", 1:length(propRep.logit))
propRep.logit <- lapply(names(propRep.logit), function(name){
  x <- propRep.logit[[name]]
  x$Replicates <- name
  return(x)
})
prop.trans.logit <- do.call(rbind, propRep.logit)
res.logit <- propellerReplicatesBlock(prop.trans.logit)

dir_path = paste0(dir, "/propeller_limma_replicates_block/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "logit_", set_level, ".rds")
saveRDS(res.logit, save_file_name)

propRep.asin <- lapply(thetaRep, function(x){
  TransformReplicates(x, "asin")
})
names(propRep.asin) <- paste0("Replicates", 1:length(propRep.asin))
propRep.asin <- lapply(names(propRep.asin), function(name){
  x <- propRep.asin[[name]]
  x$Replicates <- name
  return(x)
})
prop.trans.asin <- do.call(rbind, propRep.asin)
res.asin <- propellerReplicatesBlock(prop.trans.asin)

save_file_name <- paste0(dir_path, "asin_", set_level, ".rds")
saveRDS(res.asin, save_file_name)
