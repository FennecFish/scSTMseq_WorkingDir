setwd("/proj/milovelab/wu/scLDAseq")
# library(splatter)
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
library(sva)

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

# # batch effect removal using combat seq
# adjusted_counts <- ComBat_seq(counts(sims), batch=sims$Batch, group=NULL)
# counts(sims) <- adjusted_counts
# cat("Batch Effect Removed \n")

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
# scSTM.mod <- stm(sce = sims, K = ngroup, prevalence = ~time, init.type = "Spectral")
scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = NULL, content = NULL,
                         N = 10, ts_runs = 50, random_run = 50,
                         max.em.its = 100, net.max.em.its = 15) #, sample = "Batch",)

all_values <- unlist(scSTM.mod$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(res, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_combat_f_nc_null/", 
                           "scSTM_", set_level, ".rds"))
# args <- prepsce(sims)
# # cat("Preparing data...\n")
# documents <- args$documents
# vocab <- args$vocab
# data <- args$meta
# 
# r.file <- paste0("../stm/R/",list.files("../stm/R/"))
# sapply(r.file, source)
# sourceCpp("../stm/src/STMCfuns.cpp")
# 
# test <- stm(documents = documents, vocab = vocab, data = data,
#             K = ngroup, prevalence = ~time, init.type = "Spectral")

