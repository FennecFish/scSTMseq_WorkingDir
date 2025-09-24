# Goal: This script is to generate a single patient simulation
# both pre and post timepoints proportion is generated from a Logistic Normal
# Gamma is drawn from a multivariate normal
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library("scater")
library(SingleCellExperiment)
library(MASS)
library(VariantAnnotation)
library(checkmate)
library(MCMCpack)

args <- commandArgs(trailingOnly = TRUE)
sim_index <- as.numeric(args[1])
seed <- as.integer(Sys.time() + sim_index * runif(1, 2,10))
set.seed(seed)

# Define parameters
batch.rmEffect = TRUE
cancerCellGroup = 2
# numCellType <- c(5, 10, 15)
numCellType <- c(10,12)
# gamma_sd <- c(0.1, 0.15, 0.3, 0.4, 0.8)
# gamma_sd <- c(0, 0.3, 0.5, 0.8, 1, 1.5)
gamma_sd <- 0
# EffectSize <- c(-0.5, -0.8, -1)
nSample <- 10
nTimepoints <- 2
A = 1 # number of covariates, only 1 to show timeEffect
nGenes = 3000
de.prob <- 0.3
de.facLoc <- 0.5
batchCells <- rep(c(250,250), each = nSample)
batch.facLoc <- runif(nSample, min = 0, max = 0.5)
batch.facLoc <- rep(batch.facLoc, times =2)

################################################################################
############################## Useful Functions ################################
################################################################################
generate_theta <- function(nSample, nTimepoints = 2, nCellType, mean, sd, Response = 1){
  Simplex <- nCellType - 1
  # sample_ids <- rep(1:nSample, each = nTimepoints) 
  Timepoint <- rep(c(0, 1), nSample)
  
  # Create design matrix for each sample/time combination
  X <- data.frame(Timepoint = Timepoint)
  if(nSample > 1){
    sample_ids <- rep(1:nSample, each = nTimepoints) 
    rownames(X) <- paste(paste0("Sample", sample_ids),
                         ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  }else{
    rownames(X) <- paste(ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  }
  
  # this maps simplex to nCEllType
  map_to_simplx <- function(x) {
    exp(x - log(sum(exp(x))))
  }

  # # randomly generate a base simplex
  # simplex_t1 <- rnorm(simplex)
  # index <-  which(simplex_t1 == max(simplex_t1))
  # # t2 will increase or decrease on the first cell type by effect size
  # simplex_t2 <- simplex_t1
  # simplex_t2[index] <- simplex_t2[index]*(1 + EffectSize)
  # 
  # mu <- data.frame()
  # if(nSample > 1){
  #   # simulate patient random effect
  #   psi <- mvrnorm(n = nSample, mu = rep(0, Simplex),  Sigma= diag(rep(1, Simplex)))
  #   
  #   # adding random effect eta = mu + psi
  #   for (i in 1:nSample) {
  #     temp <- rbind(simplex_t1 + psi[i, ], simplex_t2 + psi[i, ])
  #     mu <- rbind(mu, temp)
  #   }
  #   param <- expand.grid(paste0("Sample", 1:nSample), c("t1", "t2")) %>%
  #     arrange(Var1, Var2)
  #   rownames(mu) <- apply(param, 1, function(x) paste(x, collapse = "_"))
  # }
  
  # Generate Gamma from a Normal Distribution
  # Using normal cause we only have one covariate
  gamma <- vector(mode = "list")
  for (i in 1:Simplex) {
    gamma[[i]] <- rnorm(n = 1, mean = mean[i], sd = sd[[i]])
  }
  gamma <- do.call(cbind, gamma)
  colnames(gamma) <- c(paste0("K", 1:ncol(gamma)))
  rownames(gamma) <- c("Timepoint")

  mu <- t(t(as.matrix(gamma))  %*% t(as.matrix(X)))

  if(nSample > 1){
    # simulate patient random effect
    psi <- mvrnorm(n = nSample, mu = rep(0, Simplex),  Sigma= diag(rep(1, Simplex)))

    # adding random effect eta = mu + psi
    for (i in 1:nSample) {
      mu[paste0("Sample",i,"_t1"), ] <- mu[paste0("Sample",i,"_t1"), ] + psi[i, ]
      mu[paste0("Sample",i,"_t2"), ] <- mu[paste0("Sample",i,"_t2"), ] + psi[i, ]
    }
  }
  # transfer the proportion into simplex.
  # adding 0 to every sample/time for identifiability in the model
  eta <- cbind(mu, 0)
  colnames(eta) <- paste0("K", 1:ncol(eta))
  
  theta <- t(apply(eta, 1, map_to_simplx))
  theta <- list(t1 = theta[grep("t1", rownames(theta)), ], t2 = theta[grep("t2", rownames(theta)), ])
  
  # return(list(theta = theta, index = index, psi = psi))

  if(nSample > 1){
    return(list(theta = theta, gamma = gamma, psi = psi))
  }else{
    return(list(theta = theta, gamma = gamma))
  }
}

# take the same baseline_proportion but varying effectSize to simulate data
delta_sim <- function(true_param, seed = seed, nSample = nSample, 
                      nGenes = nGenes, nCellType = nCellType, 
                      de.prob = de.prob, de.facLoc = de.facLoc, 
                      batchCells = batchCells, batch.facLoc = batch.facLoc,
                      cancerCellGroup = NULL,
                      batch.rmEffect = TRUE){
  params <- newSplatParams()
  params <- setParams(params, nGenes = nGenes,
                      group.prob = true_param$theta,
                      de.prob = de.prob, de.facLoc = de.facLoc,
                      batchCells=batchCells, batch.facLoc = batch.facLoc,
                      seed = seed)
  sims <- splatSimulate(params, method = "groups",
                        verbose = F, batch.rmEffect = batch.rmEffect,
                        cancerCellGroup = cancerCellGroup, 
                        true_param = true_param)
  return(sims)
}

# Write a function to check sims proportion change
proportion_check <- function(sims){
  sims_df <- as.data.frame(colData(sims))
  
  # Filter for Time1 and Time2
  sims_time1 <- sims_df[sims_df$Time == "Time1", ]
  sims_time2 <- sims_df[sims_df$Time == "Time2", ]
  
  # Count the number of cells in each group for Time1 and Time2
  group_counts_time1 <- table(sims_time1$Sample, sims_time1$Group)
  group_counts_time2 <- table(sims_time2$Sample, sims_time2$Group)
  
  # Convert the counts to proportions (relative frequency)
  prop_time1 <- prop.table(group_counts_time1, margin = 1)
  prop_time2 <- prop.table(group_counts_time2, margin = 1)
  
  # Calculate the proportion change for each group between Time1 and Time2
  proportion_change <- (prop_time2 - prop_time1)/prop_time1
  
  return(list(t1 = prop_time1, t2= prop_time2))
  # Convert to data frame for better readability
  # proportion_change_df <- as.data.frame(as.table(proportion_change))
}


################################################################################
############################### Simulations ####################################
################################################################################

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
if(batch.rmEffect){save_batch <- "noBatch"} else{save_batch <- "Batch"}
if(is.null(cancerCellGroup)){save_cancer <- "StromalCell"} else{save_cancer <- "CancerCell"}

r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

param_dat <- expand.grid(numCellType, gamma_sd)
colnames(param_dat) <- c("numCellType", "gamma_sd")

# #############################################################
for (i in 1:nrow(param_dat)){
  nCellType <- param_dat[i,1]
  gamma_sd_tmp <- param_dat[i,2]
  
  simplex <- nCellType- 1 
  mean = matrix(rep(0, simplex), nrow = 1)
  sd <- replicate(simplex, diag(gamma_sd_tmp, A), simplify = FALSE)
  if(gamma_sd_tmp == 0){
    type = "NullModel"
  }else{
    type = paste0("HighVar", gamma_sd_tmp)
  }

  true_param <- generate_theta(nSample = nSample, nTimepoints = 2, nCellType = nCellType, mean = mean, sd = sd)
  sims <- delta_sim(true_param = true_param, seed = seed, nSample = nSample,
                    nGenes = nGenes, nCellType = nCellType,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    batchCells = batchCells, batch.facLoc = batch.facLoc,
                    cancerCellGroup = cancerCellGroup,
                    batch.rmEffect = batch.rmEffect)
  
  # save_response <- paste(round(ResponseEffect, 1), collapse = ".")
  # save_time <- paste(round(TimeEffect, 1), collapse = ".")
  dir_path <- paste0(dir, "nSample", nSample,
                     "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  save_file_name <- paste0(dir_path, "sims_", seed, "_", type, ".rds")
  saveRDS(sims, file = save_file_name)
  # print(proportion_check(sims))
  # print((true_param$theta$t2-true_param$theta$t1)/true_param$theta$t1)
  rm(sims)
  cat("Cell Type is", param_dat[i,1], "and the sd is", param_dat[i,2], "\n")
  # sims <- sims[,sims$Sample == c("Sample1", "Sample2")]
  # sims <- logNormCounts(sims)
  # sims <- runPCA(sims)
  # plotPCA(sims, colour_by = "Batch", shape_by = "Sample")
  # plotPCA(sims, colour_by = "Batch", shape_by = "Group")
  # plotPCA(sims, colour_by = "Sample")
  
  # #############################################################
  # mean = matrix(rep(0, simplex*A), nrow = A)
  # sd <- replicate(simplex, diag(gamma_sd_tmp, A), simplify = FALSE)
  # type = "HighVar"
  # 
  # true_param <- generate_theta(nSample = nSample, nTimepoints = 2, nCellType = nCellType, mean = mean, sd = sd)
  # sims <- delta_sim(true_param = true_param, seed = seed, nSample = nSample,
  #                   nGenes = nGenes, nCellType = nCellType,
  #                   de.prob = de.prob, de.facLoc = de.facLoc,
  #                   batchCells = batchCells, batch.facLoc = batch.facLoc,
  #                   cancerCellGroup = cancerCellGroup,
  #                   batch.rmEffect = batch.rmEffect)
  # 
  # # save_response <- paste(round(ResponseEffect, 1), collapse = ".")
  # # save_time <- paste(round(TimeEffect, 1), collapse = ".")
  # 
  # saveRDS(sims, file = paste0(dir, "nSample", nSample,
  #                             "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/",
  #                             "sims_", seed, "_", type, ".rds"))
  # print(proportion_check(sims))
  # print((true_param$theta$t2-true_param$theta$t1)/true_param$theta$t1)
  # rm(sims)
  
}

