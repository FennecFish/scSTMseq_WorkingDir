# This script is for sensitivity analysis
# assessing how many iterations are needed 
setwd("/proj/milovelab/wu/scLDAseq")
library(fastTopics)
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
library(mclust)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
file_name <- args[2]
design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
cat("Design is", design, "\n")
file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat(file_name, "\n")
gc <- as.integer(args[3])
iter <- as.integer(args[4])
num.init <- as.integer(args[5])
cat("Number of Iteration is", iter, "with num.init", num.init,"\n")
set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
sims <- readRDS(paste0(dir, "/sims/", file_name))
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

# test 
# sims <- sims[1:100,2400:2600]

# nsample <- length(unique(sims$Sample))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

# test <- selectModel(sce = sims, sample = "Sample",
#                      K = ngroup, prevalence = ~Time, content = ~Sample,
#                      gamma.prior = "Pooled",
#                      N = 5, ts_runs = 2, random_run = 2,
#                      max.em.its = 100, net.max.em.its = 2)

N <- 2*num.init + 1
scSTM.mod <- selectModel_parallel(sce = sims, sample = "Sample",
                                  K = ngroup, prevalence = ~Time, content = NULL,
                                  gamma.prior = "Pooled",
                                  N = N, ts_runs = num.init, random_run = num.init,
                                  max.em.its = 100, net.max.em.its = iter, gc = gc)

# msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
dir_path <- paste0(dir, "/scSTM_Sensitivity_Pooled_Content_Prevalence_Time/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "scSTM_iter", iter, "_init", num.init, "_", set_level, "_ARI.rds")
saveRDS(scSTM.mod, file = save_file_name)
