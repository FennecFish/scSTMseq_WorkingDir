# This script is designed to compare results from stm and scLDAseq
# when running on individual patient from Anti-PD1 data

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

patient <- read.csv("doc/patientID.csv")
patientID <- patient$x[8]
sims.sub <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/PD1-data/", patientID, ".rds"))

# stm estimate
library(stm)
stm_res <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/res/PD1-res/stm_", patientID, ".rds"))
K <- stm_res$settings$dim$K
meta <- colData(sims.sub) %>% data.frame()
stm_est <- estimateEffect(1:K ~ timepoint, stm_res, meta=meta, 
                       uncertainty="Global")
summary(stm_est)
plot(stm_est, "timepoint", xlim = c(-0.2,0.5),
     main = paste0("Estimated Effect of Time as a Covariate to Cell-Topic Proportion for STM in Patient", patientID),
     method="difference",cov.value1="Pre",cov.value2="On",
     model = stm_res,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")

#scLDAseq estimate
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

sc_res <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/res/PD1-res/scLDAseq_", patientID, ".rds"))
K <- sc_res$settings$dim$K
sc_est <- estimateEffect(1:K ~ timepoint, 
                               stmobj = sc_res, 
                               meta= meta, uncertainty = "Global")
summary(sc_est)
plot(sc_est, "timepoint", xlim = c(-0.2,0.5),
     main = paste0("Estimated Effect of Time as a Covariate to Cell-Topic Proportion for scLDAseq in Patient ", patientID),
     method="difference",cov.value1="Pre",cov.value2="On",
     model = sc_res,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")
