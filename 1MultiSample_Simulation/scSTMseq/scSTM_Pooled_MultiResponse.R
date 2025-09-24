# this script is to run scSTMseq with Time and Response VAriable
# with content
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

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
file_name <- args[2]
design <- sub(".*/MultiResponse/([^/]+)/.*", "\\1", file_name)
cat("Design is", design, "\n")
file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat(file_name, "\n")
gc <- as.integer(args[3])

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

# nsample <- length(unique(sims$Sample))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

# selectModel(sce = sims, sample = NULL,
#                      K = ngroup, prevalence = ~Time, content = NULL,
#                      gamma.prior = "Pooled",
#                      N = 1, ts_runs = 1, random_run = 1,
#                      max.em.its = 1, net.max.em.its = 2)

scSTM.mod <- selectModel_parallel(sce = sims, sample = "Sample",
                                  K = ngroup, prevalence = ~Time*Response, content = ~Sample,
                                  gamma.prior = "Pooled",
                                  N = 5, ts_runs = 30, random_run = 30,
                                  max.em.its = 100, net.max.em.its = 15, gc = gc)

# msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
dir_path <- paste0(dir, "/scSTM_Pooled_Content_Sample_Prevalence_Time/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "scSTM_", set_level, ".rds")
saveRDS(scSTM.mod, file = save_file_name)
