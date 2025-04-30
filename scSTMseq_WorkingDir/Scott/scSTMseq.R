setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
# library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(dplyr)
library(cluster)
library(scater)
library(scran)

# read in the simulated data
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

# extract simulation setup from the name
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)

# read in teh data
sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/composition_change/", file_name))

# quick qc
sims <- quickPerCellQC(sims, filter=TRUE)
### remove genes with count 0 
sims <- sims[rowSums(counts(sims)) != 0,]
# #### feature selection #####
sims <- scuttle::logNormCounts(sims)
dec.p2 <- modelGeneVar(sims)
# feature selection
p2.chosen <- getTopHVGs(dec.p2, n=2000)
sims <- sims[p2.chosen,]

nsample <- length(unique(sims$Batch))
ngroup <- length(unique(sims$Group))

# read in all the R files and Rcpp files in the R folder
r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()

# this is used to perform scSTMseq directly without content model
res <- multi_stm(sce = sims,
                 K = ngroup, prevalence = ~time, content = NULL,
                 sample = "Batch",
                 init.type= "Random",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000),
                 emtol=1e-5)

# with content model. This will influence the derivation of beta
res <- multi_stm(sce = sims,
                 K = ngroup, prevalence = ~time, content = ~Batch,
                 sample = "Batch",
                 init.type= "Random",
                 gamma.prior= "Pooled",
                 kappa.prior= "L1",
                 control = list(gamma.maxits=3000),
                 emtol=1e-5)


# this is used to select best performance model with content model
scSTM.mod <- selectModel(sce = sims,
                         K = ngroup, prevalence = ~time, content = ~Batch,
                         sample = "Batch", N = 5, ts_runs = 20, random_run = 20)

all_values <- unlist(scSTM.mod$bound)
max_value <- max(all_values)
max_position_in_vector <- which(all_values == max_value)
res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scSTMseq (%d seconds). \n", floor((proc.time()-t1)[3]))
saveRDS(res, file = paste0("scSTM_", set_level, ".rds"))
