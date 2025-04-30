# This script is designed to run stm on Anti-PD1 data
setwd("/proj/milovelab/wu/scLDAseq/scLDAseq")
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


dat <- readRDS("/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/raw/cohort1_filtered.rds")
# subdat <- subset(dat, subset = expansion=="E")
subdat <- subset(dat, subset = patient_id =="BIOKEY_26") #28,3 31
sims <- as.SingleCellExperiment(subdat)
rm(dat)
rm(subdat)
cat("Completed Transition to sce, starting QC \n")

#### QC ######
sims <- quickPerCellQC(sims,filter=TRUE)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2000)
sims.sub <- sims[p2.chosen,]
rm(sims)
# saveRDS(sims.sub, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_E_only_2000genes.rds")
###########################################################
################## scLDAseq ###############################
###########################################################
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

stm_dat <- prepsce(sims.sub)

# saveRDS(stm_dat, file = "/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/stm_prepsce_anti-PD1_E_only_2000genes.rds")
cat("Prepreation Completed")

K <- 40

# stm_dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/stm_prepsce_anti-PD1.rds")
prevalence <- as.formula(~stm_dat$meta$timepoint)
content <- NULL
sce <- stm_dat$sce
documents  <- stm_dat$documents
vocab <- stm_dat$vocab
data <- stm_dat$meta
sample <- "patient_id"

res.stm <- multi_stm(documents = documents, vocab = vocab,
                     K = K, prevalence = prevalence, content = NULL,
                     data = data,
                     sce = sce,
                     sample = sample,
                     init.type= "Spectral",
                     gamma.prior= "Pooled",
                     kappa.prior= "L1",
                     control = list(gamma.maxits=3000))

saveRDS(res.stm, file = "res/Anti-PD1_BIOKEY15_1000genes_K40_stmRes.rds")
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))