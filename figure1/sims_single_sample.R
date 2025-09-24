# this script is to generate simulated data for figure 1
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(splatter)
library("scater")
library(SingleCellExperiment)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))

##################################### level 1 #####################################
# 2000 genes, 2500 cells, 5 cell types, 
# each with differential probability of 0.5, and de.facLoc = 2
nsample <- 1
nCellType <- 5
# batchCells <- 5000
batchCells <- rep(1500, nsample)
nGenes <- 2000
dropout.type = "batch"
dropout.mid <- 0.05
dropout.shape <- -1
de.prob <- runif(nCellType, min = 0.3, max = 0.6)
de.facLoc <- runif(nCellType, 1.5, 2.5)

comp <- c(0.2, 0.2, 0.2, 0.2, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s1 <- sims
# sim3 <- logNormCounts(s1)
# sim3 <- runPCA(sim3)
# plotPCA(sim3, colour_by = "Group")

### negative control
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s2 <- sims
# sim4 <- logNormCounts(s2)
# sim4 <- runPCA(sim4)
# plotPCA(sim4, colour_by = "Group")

s1$time <- 1
s1$Cell <-paste0(s1$Cell, "_1")
rownames(colData(s1)) <- s1$Cell
s2$time <- 2
s2$Cell <-paste0(s2$Cell, "_2")
rownames(colData(s2)) <- s2$Cell

s <- cbind(s1, s2)
# s <- logNormCounts(s)
# s <- runPCA(s)
# plotPCA(s, colour_by = "time")

cat("Generation Completed for L1.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                            "sims_", seed, "_neg_L1_c5.rds"))

#### positive change
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

comp <- c(0.35, 0.1, 0.1, 0.25, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s3 <- sims

s3$time <- 2
s3$Cell <-paste0(s3$Cell, "_2")
rownames(colData(s3)) <- s3$Cell

s <- cbind(s1, s3)
# s <- logNormCounts(s)
# s <- runPCA(s)
# plotPCA(s, colour_by = "time")

# df <- colData(s) %>% 
#   as.data.frame() %>% 
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0))
# total_1 <- sum(df$`1`)
# total_2 <- sum(df$`2`)
# result <- df %>%
#   mutate(
#     Proportion_1 = `1` / total_1,
#     Proportion_2 = `2` / total_2,
#     Ratio = Proportion_1 / Proportion_2
#   )

cat("Generation Completed for L1-pos.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                         "sims_", seed, "_pos_L1_c5.rds"))


##################################### level 2 #####################################
# 2000 genes, 2500 cells, 5 cell types, 
# each with differential probability of 0.3, and de.facLoc = 1
nsample <- 1
nCellType <- 5
batchCells <- rep(1500, nsample)
nGenes <- 2000
dropout.type = "batch"
dropout.mid <- 0.05
dropout.shape <- -1
de.prob <- runif(nCellType, min = 0.1, max = 0.3)
de.facLoc <- runif(nCellType, 0.8, 1.5)

comp <- c(0.2, 0.2, 0.2, 0.2, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s1 <- sims
# sim3 <- logNormCounts(s1)
# sim3 <- runPCA(sim3)
# plotPCA(sim3, colour_by = "Group")

### negative control #####
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s2 <- sims
# sim4 <- logNormCounts(s2)
# sim4 <- runPCA(sim4)
# plotPCA(sim4, colour_by = "Group")

s1$time <- 1
s1$Cell <-paste0(s1$Cell, "_1")
rownames(colData(s1)) <- s1$Cell
s2$time <- 2
s2$Cell <-paste0(s2$Cell, "_2")
rownames(colData(s2)) <- s2$Cell

s <- cbind(s1, s2)
# s <- logNormCounts(s)
# s <- runPCA(s)
# plotPCA(s, colour_by = "Group")

cat("Generation Completed for L2-neg.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                         "sims_", seed, "_neg_L2_c5.rds"))

#### positive change #####
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

comp <- c(0.25, 0.1, 0.1, 0.25, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s3 <- sims
# sim1 <- logNormCounts(s3)
# sim1 <- runPCA(sim1)
# plotPCA(sim1, colour_by = "Group")

s3$time <- 2
s3$Cell <-paste0(s3$Cell, "_2")
rownames(colData(s3)) <- s3$Cell

s <- cbind(s1, s3)
# s <- logNormCounts(s)
# s <- runPCA(s)
# plotPCA(s, colour_by = "Group")

# df <- colData(s) %>% 
#   as.data.frame() %>% 
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0))
# total_1 <- sum(df$`1`)
# total_2 <- sum(df$`2`)
# result <- df %>%
#   mutate(
#     Proportion_1 = `1` / total_1,
#     Proportion_2 = `2` / total_2,
#     Ratio = Proportion_1 / Proportion_2
#   )

cat("Generation Completed for L2-pos.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                         "sims_", seed, "_pos_L2_c5.rds"))

##################################### level 3 #####################################
# 2000 genes, 2500 cells, 5 cell types, 
# each with differential probability of 0.2, and de.facLoc = 0.05
nsample <- 1
nCellType <- 5
# batchCells <- 5000
batchCells <- rep(1500, nsample)
nGenes <- 2000
dropout.type = "batch"
dropout.mid <- 0.05
dropout.shape <- -1
de.prob <- runif(nCellType, min = 0.05, max = 0.2)
de.facLoc <- runif(nCellType, 0.05, 0.2)

comp <- c(0.2, 0.2, 0.2, 0.2, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- newSplatParams()
params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s1 <- sims
# sim3 <- logNormCounts(s1)
# sim3 <- runPCA(sim3)
# plotPCA(sim3, colour_by = "Group")

### negative control #####
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s2 <- sims
# sim4 <- logNormCounts(s2)
# sim4 <- runPCA(sim4)
# plotPCA(sim4, colour_by = "Group")

s1$time <- 1
s1$Cell <-paste0(s1$Cell, "_1")
rownames(colData(s1)) <- s1$Cell
s2$time <- 2
s2$Cell <-paste0(s2$Cell, "_2")
rownames(colData(s2)) <- s2$Cell

s <- cbind(s1, s2)
# s <- logNormCounts(s)
# s <- runPCA(s)
# plotPCA(s, colour_by = "Group")

cat("Generation Completed for L3-neg.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                         "sims_", seed, "_neg_L3_c5.rds"))

#### positive change #####
params <- newSplatParams()
bcv.common <- 0.2
lib.loc <- 13
lib.scale <- 0.25

comp <- c(0.35, 0.1, 0.1, 0.25, 0.2) * sum(batchCells)
de.comp <- mvrnorm(n = nsample, mu = comp, Sigma = diag(x=200, nrow =5))
group.prob <- de.comp/sum(de.comp)

params <- setParams(params, group.prob = group.prob,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    nGenes = nGenes, batchCells=batchCells,
                    seed = seed,
                    bcv.common = bcv.common,
                    lib.loc = lib.loc, lib.scale = lib.scale)
sims <- splatSimulate(params, method = "groups",
                      verbose = TRUE, batch.rmEffect = FALSE)
s3 <- sims

s3$time <- 2
s3$Cell <-paste0(s3$Cell, "_2")
rownames(colData(s3)) <- s3$Cell

s <- cbind(s1, s3)
# sim6 <- logNormCounts(s)
# sim6 <- runPCA(sim6)
# plotPCA(sim6, colour_by = "Group")

# df <- colData(s) %>% 
#   as.data.frame() %>% 
#   count(Batch, Group, time) %>%
#   pivot_wider(names_from = time, values_from = n, values_fill = list(n = 0))
# total_1 <- sum(df$`1`)
# total_2 <- sum(df$`2`)
# result <- df %>%
#   mutate(
#     Proportion_1 = `1` / total_1,
#     Proportion_2 = `2` / total_2,
#     Ratio = Proportion_1 / Proportion_2
#   )

cat("Generation Completed for L3-pos.\n")

saveRDS(s, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V2_single/sims/", 
                         "sims_", seed, "_pos_L3_c5.rds"))