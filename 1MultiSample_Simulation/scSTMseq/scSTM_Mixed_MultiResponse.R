# this script is to run scSTMseq with filtered gene
# without content
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stm)
library(scater)
library(scran)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")
gc <- as.integer(args[2])

set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_MultiResponse/nSample6_nCellType5_noBatch_StromalCell/"
sims <- readRDS(paste0(dir, "sims/", file_name))
# sims$Batch_ID <- paste0(sims$Sample, "_", sims$Time)
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1500)
sims <- sims[p2.chosen,]

# nsample <- length(unique(sims$Sample))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
# scSTM.mod <- stm(sce = sims, K = ngroup, prevalence = ~time, init.type = "Spectral")
# scSTM.mod <- selectModel(sce = sims, sample = NULL,
#                          K = ngroup, prevalence = ~Time*Response + (1|Sample), content = NULL,
#                          gamma.prior = "LinearMixed",
#                          N = 1, ts_runs = 1, random_run = 1,
#                          max.em.its = 5, net.max.em.its = 2)
# scSTM.mod <- selectModel(sce = sims, sample = NULL,
#                          K = ngroup, prevalence = ~Time, content = NULL,
#                          gamma.prior = "LinearRegression",
#                          N = 1, ts_runs = 1, random_run = 1,
#                          max.em.its = 2, net.max.em.its = 15)


# scSTM.mod <- selectModel_parallel(sce = sims, sample = "Sample",
#                                   K = ngroup, prevalence = ~Time*Response, content = ~Sample,
#                                   N = 5, ts_runs = 30, random_run = 30,
#                                   max.em.its = 100, net.max.em.its = 15, gc = gc)

scSTM.mod <- selectModel_parallel(sce = sims, sample = NULL,
                                  K = ngroup, prevalence = ~Time*Response + (1|Sample), content = NULL,
                                  gamma.prior = "LinearMixed",
                                  N = 5, ts_runs = 30, random_run = 30,
                                  max.em.its = 100, net.max.em.its = 15, gc = gc)

msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(scSTM.mod, file = paste0(dir, "scSTM_LinearMixed_noContent_Prevalence_TimeandResponse/", 
                                 "scSTM_", set_level, ".rds"))

