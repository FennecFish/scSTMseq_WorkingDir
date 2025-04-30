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

sims_scRNA <- function(nsample = 1, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed){
  
  batchCells <- rep(200, nsample)
  nGenes <- 2000
  dropout.type = "batch"
  dropout.mid <- 0.05
  dropout.shape <- -1
  sd <- batchCells[1] * 0.05
  # comp <- c(0.2, 0.2, 0.2, 0.2, 0.2) * sum(batchCells)
  comp <- rep(1/nCellType, nCellType) * sum(batchCells)
  de.comp <- mvrnorm(n = 1, mu = comp, Sigma = diag(x=sd, nrow =nCellType))
  group.prob <- de.comp/sum(de.comp)
  
  params <- newSplatParams()
  params <- setParams(params, group.prob = group.prob,
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      nGenes = nGenes, seed = seed,
                      batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale)
  
  sims <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = FALSE)
  s1 <- sims
  # sim1 <- logNormCounts(s1)
  # sim1 <- runPCA(sim1)
  # plotPCA(sim1, colour_by = "Batch")
  
  ### negative control
  de.comp <- mvrnorm(n = 1, mu = comp, Sigma = diag(x=sd, nrow =nCellType))
  group.prob <- de.comp/sum(de.comp)
  
  params <- newSplatParams()
  params <- setParams(params, group.prob = group.prob,
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      nGenes = nGenes, seed = seed,
                      batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale)
  sims <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = FALSE)
  s2 <- sims
  # sim2 <- logNormCounts(s2)
  # sim2 <- runPCA(sim2)
  # plotPCA(sim2, colour_by = "Group")
  
  s1$time <- 1
  s1$Cell <-paste0(s1$Cell, "_1")
  rownames(colData(s1)) <- s1$Cell
  s2$time <- 2
  s2$Cell <-paste0(s2$Cell, "_2")
  rownames(colData(s2)) <- s2$Cell
  
  sim_neg_dat <- cbind(s1, s2)
  # sim_ng <- logNormCounts(sim_neg_dat)
  # sim_ng <- runPCA(sim_ng)
  # plotPCA(sim_ng, colour_by = "Batch")
  # plotPCA(sim_ng, colour_by = "time")
  
  # 
  #### positive change
  # if(nCellType == 5){change.comp <- rep(1/nCellType, nCellType) + c(0.2, -0.15, -0.1, 0.05, 0.2)} 
  # if(nCellType == 9){change.comp <- rep(1/nCellType, nCellType) + c(0.05, -0.05, 0.08, -0.03, -0.05, 0.02, -0.02, 0, 0)} 
  # if(nCellType == 13){change.comp <- rep(1/nCellType, nCellType) + c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
  #                                                                    -0.04, 0)} 
  if(nCellType == 6){change.comp <- rep(1/nCellType, nCellType) + c(0.1, -0.12, -0.1, 0.05, 0.07, 0)} 
  if(nCellType == 12){change.comp <- rep(1/nCellType, nCellType) + c(0.03, -0.01, -0.04, 0.02, 0.02, -0.05, -0.01, 0.1, -0.02, -0.03, -0.01, 0)} 
  # if(nCellType == 13){change.comp <- rep(1/nCellType, nCellType) + c(0.03, -0.01, -0.02, 0.02, 0.02, -0.03, -0.01, 0.1, -0.02, -0.03, -0.01,
  #                                                                    -0.04, 0)} 
  comp <- change.comp * sum(batchCells)
  
  de.comp <- mvrnorm(n = 1, mu = comp, Sigma = diag(x=sd, nrow = nCellType))
  group.prob <- de.comp/sum(de.comp)
  
  params <- newSplatParams()
  params <- setParams(params, group.prob = group.prob,
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      nGenes = nGenes, seed = seed,
                      batchCells=batchCells, batch.facLoc = batch.facLoc, batch.facScale = batch.facScale)
  sims <- splatSimulate(params, method = "groups",
                        verbose = TRUE, batch.rmEffect = FALSE)
  s3 <- sims
  # sim3 <- logNormCounts(s3)
  # sim3 <- runPCA(sim3)
  # plotPCA(sim3, colour_by = "Group")
  
  s3$time <- 2
  s3$Cell <-paste0(s3$Cell, "_2")
  rownames(colData(s3)) <- s3$Cell
  
  sim_pos_dat <- cbind(s1, s3)
  # sim_pos <- logNormCounts(sim_pos_dat)
  # sim_pos <- runPCA(sim_pos)
  # plotPCA(sim_pos, colour_by = "time")
  # plotPCA(sim_pos, colour_by = "Batch")
  
  return(list(sim_neg_dat, sim_pos_dat))
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
  
 
}

# nsample <- c(3, 5, 10)
# nCellType <- c(5, 9, 13)

nsample <- c(5, 10)
nCellType <- c(6, 12)


nparam <- expand.grid(nsample, nCellType)
colnames(nparam) <- c("nsample", "nCellType")

for(i in 1:dim(nparam)[1]){

  nsample <- nparam[i,1]
  nCellType <- nparam[i,2]
  
  #### level 1 ####
  de.prob <- runif(nCellType, min = 0.5, max = 1)
  de.facLoc <- runif(nCellType, 2, 2.5)
  
  batch.facLoc = runif(nsample, min = 0, max = 0.1)
  batch.facScale = runif(nsample, min = 0, max = 0.1)
  
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                           "sims_", seed, "_neg_L1_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L1-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L1_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L1-pos.\n")
  rm(sims)
  
  #### level 2 ####
  de.prob <- runif(nCellType, min = 0.3, max = 0.6)
  de.facLoc <- runif(nCellType, 1, 1.5)
  
  batch.facLoc = runif(nsample, min = 0, max = 0.1)
  batch.facScale = runif(nsample, min = 0, max = 0.1)
  
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_neg_L2_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L2-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L2_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L2-pos.\n")
  rm(sims)
  #### level 3 ####
  de.prob <- runif(nCellType, min = 0.1, max = 0.3)
  de.facLoc <- runif(nCellType, 0.5, 1)
  
  batch.facLoc = runif(nsample, min = 0, max = 0.1)
  batch.facScale = runif(nsample, min = 0, max = 0.1)
  
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_neg_L3_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L3-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L3_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L3-pos.\n")
  rm(sims)
  #### level 4 ####
  de.prob <- runif(nCellType, min = 0.1, max = 0.3)
  de.facLoc <- runif(nCellType, 0.5, 1)

  batch.facLoc = runif(nsample, min = 0.1, max = 0.3)
  batch.facScale = runif(nsample, min = 0.1, max = 0.3)

  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)

  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_neg_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4-pos.\n")
  rm(sims)
  #### level 5 ####
  de.prob <- runif(nCellType, min = 0.1, max = 0.3)
  de.facLoc <- runif(nCellType, 0.1, 0.3)

  batch.facLoc = runif(nsample, min = 0.2, max = 0.5)
  batch.facScale = runif(nsample, min = 0.2, max = 0.5)

  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)

  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_neg_L5_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L5-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L5_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L5-pos.\n")
  rm(sims)
  #### level 5 ####
  de.prob <- runif(nCellType, min = 0.05, max = 0.15)
  de.facLoc <- runif(nCellType, 0.05, 0.15)
  
  batch.facLoc = runif(nsample, min = 0.2, max = 0.5)
  batch.facScale = runif(nsample, min = 0.2, max = 0.5)
  
  sims <- sims_scRNA(nsample, nCellType, de.prob, de.facLoc, batch.facLoc, batch.facScale, seed)
  
  saveRDS(sims[[1]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_neg_L4_c", nCellType, "_nsample", nsample, ".rds"))
  cat("Generation Completed for L4-neg.\n")
  saveRDS(sims[[2]], file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims/",
                                   "sims_", seed, "_pos_L4_c", nCellType, "_ns", nsample, ".rds"))
  cat("Generation Completed for L4-pos.\n")
  rm(sims)
}

# sim1 <- logNormCounts(sims[[2]])
# sim1 <- runPCA(sim1)
# plotPCA(sim1, colour_by = "Group")
# plotPCA(sim1, colour_by = "Batch")
# plotPCA(sim1, colour_by = "time")
