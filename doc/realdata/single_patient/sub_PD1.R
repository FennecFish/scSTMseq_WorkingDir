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

sims.sub <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_p3_p10_2000g.rds")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

meta <- colData(sims.sub) %>% data.frame()
sc_res <- readRDS("res/stmRes_noBatch_PD1_p3_p10_2000g.rds")
K <- sc_res$settings$dim$K
sc_est <- estimateEffect(1:K ~ timepoint, 
                         stmobj = sc_res, 
                         sampleNames = "patient_id", 
                         sampleIDs = "BIOKEY_3",
                         meta= meta, uncertainty = "Global")
summary(sc_est)
plot(sc_est, "timepoint", xlim = c(-0.3,0.3),
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion for scLDAseq in Patient BIOKEY 3",
     method="difference",cov.value1="Pre",cov.value2="On",
     model = sc_res,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")

sc_est_10 <- estimateEffect(1:K ~ timepoint, 
                         stmobj = sc_res, 
                         sampleNames = "patient_id", 
                         sampleIDs = "BIOKEY_10",
                         meta= meta, uncertainty = "Global")
summary(sc_est_10)
plot(sc_est_10, "timepoint", xlim = c(-0.3,0.3),
     main = "Estimated Effect of Time as a Covariate to Cell-Topic Proportion for scLDAseq in Patient BIOKEY 10",
     method="difference",cov.value1="Pre",cov.value2="On",
     model = sc_res,
     labeltype = "frex",
     xlab = "Proportion Difference (95% Confidence Interval) ")
