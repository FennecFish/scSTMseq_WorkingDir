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
library(Seurat)
library(cluster)
library(geometry)
library(Rtsne)
library(rsvd)
set.seed(1)

sims.sub <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_p3_p10_2000g.rds")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
# 
# stm_dat <- prepsce(sims.sub)
# 
# saveRDS(stm_dat, file = "/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/prepsce_PD1_p3_p10_2000g.rds")
# cat("Prepreation Completed")

stm_dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/prepsce_PD1_p3_p10_2000g.rds")
prevalence <- as.formula(~stm_dat$meta$timepoint)
content <- NULL
sce <- stm_dat$sce
documents  <- stm_dat$documents
vocab <- stm_dat$vocab
data <- stm_dat$meta
sample <- "patient_id"

K <- 8
res.stm <- multi_stm(documents = documents, vocab = vocab,
                     K = K, prevalence = prevalence, 
                     content = ~patient_id,
                     data = data,
                     sce = sce,
                     sample = sample,
                     init.type= "Spectral",
                     gamma.prior= "Pooled",
                     kappa.prior= "L1",
                     control = list(gamma.maxits=3000))

saveRDS(res.stm, file = "res/stmRes_Batch_PD1_p3_p10_2000g.rds")
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
cat(msg)
