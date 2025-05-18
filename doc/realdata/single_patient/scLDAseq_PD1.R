# This script is designed to run scLDAseq on Anti-PD1 data, for each patient
# this step is to compare between scLDAseq and STM
# Since the results on per patient should NOT differ much
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

setwd("/proj/milovelab/wu/scLDAseq")

args <- commandArgs(trailingOnly = TRUE)
index <- as.character(args[1])
number_part <- as.numeric(gsub("[^0-9]", "", index))
patientID <- paste0("BIOKEY_",number_part)
cat(patientID, "\n")
sims.sub <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/PD1-data/", patientID, ".rds"))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

K = 8
stm_dat <- prepsce(sims.sub)
prevalence <- as.formula(~stm_dat$meta$timepoint)
content <- NULL
sce <- stm_dat$sce
documents  <- stm_dat$documents
vocab <- stm_dat$vocab
data <- stm_dat$meta
sample <- NULL

res.stm <- multi_stm(documents = documents, vocab = vocab,
                     K = K, prevalence = prevalence, content = NULL,
                     data = data,
                     sce = sce,
                     sample = sample,
                     init.type= "Spectral",
                     gamma.prior= "Pooled",
                     kappa.prior= "L1",
                     control = list(gamma.maxits=3000))

saveRDS(res.stm, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/PD1-res/scLDAseq_", patientID, ".rds"))
