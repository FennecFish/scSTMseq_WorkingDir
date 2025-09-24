# this script is to generate simulated data for figure 1
setwd("/proj/milovelab/wu/scLDAseq")
library(dplyr)
library(scuttle)
library(tidyverse)
library(splatter)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

##################################### level 1 #####################################
# Positive and Negative Controls. One batch, 3000 genes, 5000 cells, 2 cell types, 
# each with differential probability of 0.5, and de.facLoc = 0.5
nsample <- 1
nCellType <- 3
# batchCells <- 5000
batchCells <- 2000
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
out.prob <- 0 # outlier expr prob

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob,
                    nGenes = 2000, batchCells=batchCells,
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

cat("Generation Completed for L1.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L1.rds"))

##################################### level 2 #####################################
# One batch, 3000 genes, 5000 cells, 5 cell
# types, each with differential probability of 0.5, and de.facLoc = 0.5
# Including outliers, and dropout rate in addition to level 1
nsample <- 1
nCellType <- 3
batchCells <- 2000
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
out.prob <- 0.05 # outlier expr prob
out.facScale = 0.5 # count of outlier
out.facLoc = 4 # how far from the main
dropout.type = "batch" #uses the same parameters for every cell in the same batch
dropout.shape = -1

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    out.prob = out.prob, out.facScale = out.facScale, out.facLoc = out.facLoc,
                    dropout.type = dropout.type, dropout.shape = dropout.shape,
                    nGenes = 2000, batchCells=batchCells,
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

cat("Generation Completed for L2.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L2.rds"))
##################################### level 3 ###################################
# One batch, 3000 genes, 5000 cells, 
# 5 cell types, each with differential probability of 0.5, and de.facLoc = 0.5
nsample <- 1
nCellType <- 5
batchCells <- 2000
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
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
                    nGenes = 2000, batchCells=batchCells,
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

cat("Generation Completed for L3.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L3.rds"))

##################################### level 4 ###################################
#  Six batches with batch effect. 3000 genes, 3000 cells each batch, 5 cell types, 
# each with differential probability of 0.5

nsample <- 6 
nCellType <- 5
batchCells <- rep(600, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
# rn_group <- runif(nCellType)
# group.prob <- rn_group / sum(rn_group) # cell type proportion
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.5
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


cat("Generation Completed for L4.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L4.rds"))
##################################### level 5 ###################################
# On top of experiment 4, decrease de.facLoc
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- rep(0.5, nCellType) # group de gene prob
de.facLoc <- 0.1
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


cat("Generation Completed for L5.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L5.rds"))
##################################### level 6 #################################
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.1, 0.3) # group de gene prob
de.facLoc <- 0.1
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


cat("Generation Completed for L6.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L6.rds"))
##################################### level 7 #################################
# further decrease de.facLoc
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
group.prob <- rep(1/nCellType, nCellType)
de.prob <- runif(nCellType, 0.1, 0.3) # group de gene prob
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


cat("Generation Completed for L7.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L7.rds"))
##################################### level 8 #################################
nsample <- 6 
nCellType <- 5
batchCells <- rep(500, nsample)
batch.facLoc <- runif(1,0,0.5)
batch.facScale <- runif(1,0,0.5)
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


cat("Generation Completed for L8.\n")

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L8.rds"))
##################################### level 9 ##################################
nsample <- 10
nCellType <- 8
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

saveRDS(sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", 
                            "sims_", seed, "_L9.rds"))