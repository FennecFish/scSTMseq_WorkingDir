setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
# library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(dplyr)
# library(mclust)
library(cluster)
library(scater)
library(scran)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

# file_name <- "sims_1712873833_L7.rds"

set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/", file_name))

# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)
#  plotColData(sims.qc, x = "sum", y="detected") 

### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
# #### feature selection #####
sims <- scuttle::logNormCounts(sims)

dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=1000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = ~time, content = ~Batch,
                         sample = "Batch", N = 10, ts_runs = 50, random_run = 50,
                         max.em.its = 100, net.max.em.its = 20)

all_values <- unlist(scSTM.mod$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(res, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/scSTM_f_c", 
                           "scSTM_", set_level, ".rds"))
