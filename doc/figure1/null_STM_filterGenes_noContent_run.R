setwd("/proj/milovelab/wu/scLDAseq")
set.seed(1)
library(slam)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(dplyr)
library(tibble)
library(stats)
library(scuttle)
library(scran)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

#file_name <- "sims_1716762390_neg_L1_c5.rds"
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims/", file_name))
# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
#### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

#### scSTMseq#############
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

dat <- prepsce(sims)
K <- length(unique(sims$Group))

# STM for reference
library(stm)
r.file <- paste0("../stm/R/",list.files("../stm/R/"))
sapply(r.file, source)
sourceCpp("../stm/src/STMCfuns.cpp")
res.stm <- stm(documents = dat$documents, vocab = dat$vocab,
               K = K, prevalence = NULL, content = NULL,
               data = dat$meta, 
               init.type= "Spectral",
               gamma.prior= "Pooled",
               kappa.prior= "L1",
               control = list(gamma.maxits=3000))
saveRDS(res.stm, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/STM_f_nc_null/", 
                           "STM_", set_level, ".rds"))

