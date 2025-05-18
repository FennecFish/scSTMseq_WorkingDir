# this script is to generate simulated data to find the optimal K
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1]) + 1
cat("index is", sim_index, "\n")
nCellType <- as.numeric(args[2])
cat("cellType is", nCellType, "\n")
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))
cat("seed is", seed, "\n")

nsample <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(nsample,0,0.3)
batch.facScale <- runif(nsample,0,0.3)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0, 0.2) # group de gene prob
de.facLoc <- 0.01
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "experiment" #uses the same parameters for every cell in the same batch
dropout.shape = -1

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 3000,
                    batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale,
                    seed = seed)

sims <- splatSimulate(params, method = "groups",
                      verbose = FALSE, batch.rmEffect = FALSE)

# assume number of cells are equal across time
cell_count <- matrix(colSums(table(sims$Group, sims$Batch))/2)
# generate cell type proportion at time 1
all_rn_pre <- matrix(runif(nsample * nCellType), ncol = nsample)
pre_prp <- sweep(all_rn_pre, 2, colSums(all_rn_pre), FUN="/")
#calculate count for each group/batch at time 1
pre_count <- ceiling(sweep(pre_prp, 2, t(cell_count)[1, ], "*"))
rownames(pre_count) <- paste0("Group",1:length(unique(sims$Group)))
colnames(pre_count) <- paste0("Batch", 1:length(unique(sims$Batch)))

sampled_data <- colData(sims) %>%
  data.frame() %>%
  group_by(Group, Batch) %>%
  mutate(time = 2) %>%
  ungroup()

# randomly assign time = 1 based on the count in pre_count
for (i in 1:ncol(pre_count)) {
  batch_name <- colnames(pre_count)[i]
  for (j in 1:nrow(pre_count)) {
    group_name <- rownames(pre_count)[j]
    sampled_data <- sampled_data %>%
      group_by(Group, Batch) %>%
      mutate(time = ifelse(
        Group == group_name & Batch == batch_name &
          row_number() %in% sample(row_number(), min(pre_count[j,i], n())), 1, time)) %>%
      ungroup()
  }
}
# assign the new time
sims$time <- sampled_data$time


cat("Generation Completed for L9.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/optimalK/",
                            "sims_", nCellType,"_", seed, ".rds"))
