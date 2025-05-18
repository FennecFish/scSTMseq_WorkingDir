# This script is designed to run stm on Anti-PD1 data
setwd("/proj/milovelab/wu/scLDAseq")
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)
library(geometry)
library(Seurat)
library(Rtsne)
library(rsvd)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
index <- as.character(args[1])
number_part <- as.numeric(gsub("[^0-9]", "", index))
patientID <- paste0("BIOKEY_",number_part)
cat(patientID, "\n")

dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/raw/cohort1_filtered.rds")
subdat <- subset(dat, subset = patient_id == patientID)
sims <- as.SingleCellExperiment(subdat)

#### QC ######
sims <- quickPerCellQC(sims,filter=TRUE)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2000)
sims.sub <- sims[p2.chosen,]
rm(sims)
saveRDS(sims.sub, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/PD1-data/", patientID, ".rds"))
# sims.sub <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1-data/BIOKEY_10.rds")
library(stm)
x <- t(counts(sims.sub))
stm_out <- readCorpus(x,type="Matrix")
stm_out$vocab <- unique(colnames(x))
stm_out$meta <- data.frame("cellname" = unique(rownames(x)),
                           "time" = sims.sub$timepoint)
processed <- prepDocuments(stm_out$documents, stm_out$vocab, stm_out$meta)
PrevFit_nb <- stm(processed$documents, processed$vocab, K=8, 
                  prevalence=~time, 
                  data=processed$meta, init.type="Spectral", 
                  seed=8458159)

saveRDS(PrevFit_nb, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/PD1-res/stm_", patientID, ".rds"))
