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

###############################################################################################
################### the following script is for Response Patients Only ##########################
###############################################################################################
dat <- readRDS("/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/raw/cohort1_filtered.rds")
patient <- paste0("BIOKEY_",c(4,6,15,16))
subdat <- subset(dat, subset = patient_id %in% patient)
sims <- as.SingleCellExperiment(subdat)
rm(dat)
rm(subdat)
cat("Completed Transition to sce, starting QC \n")

#### QC ######
sims <- quickPerCellQC(sims)
#### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims.sub <- sims[p2.chosen,]
rm(sims)
saveRDS(sims.sub, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_sub_mix_response_samples_1000g.rds")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()



stm_dat <- prepsce(sims.sub)

saveRDS(stm_dat, file = "/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/stm_prepsce_anti-PD1_mix_1000g.rds")
cat("Prepreation Completed")

# stm_dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/stm_prepsce_anti-PD1_mix.rds")


# stm_dat <- readRDS("/work/users/e/u/euphyw/sc_cancer_proj/data/Anti-PD1/stm_prepsce_anti-PD1.rds")
prevalence <- as.formula(~stm_dat$meta$timepoint)
content <- NULL
sce <- stm_dat$sce
documents  <- stm_dat$documents
vocab <- stm_dat$vocab
data <- stm_dat$meta
sample <- "patient_id"

K <- 18
res.stm <- multi_stm(documents = documents, vocab = vocab,
                     K = K, prevalence = prevalence, content = NULL,
                     data = data,
                     sce = sce,
                     sample = sample,
                     init.type= "Spectral",
                     gamma.prior= "Pooled",
                     kappa.prior= "L1",
                     control = list(gamma.maxits=3000))

saveRDS(res.stm, file = "res/Anti-PD1_mix_K18_stmRes.rds")
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
cat(msg)
