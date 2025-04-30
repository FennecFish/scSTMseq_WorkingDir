# Goal: This script is to generate multiple patient simulation
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
batch.rmEffect = FALSE
cancerCellGroup = 2
# numCellType <- c(5, 10, 15)
# gamma_sd <- c(0, 0.1, 0.16, 0.23, 0.3, 0.4)
# nSample <- 20
numCellType <- 5
gamma_sd <- c(1, 2, 4)
nSample <- 20
nTimepoints <- 2

A = 3 
nGenes = 3000
de.prob <- 0.3
de.facLoc <- 0.5
batchCells <- rep(c(250,250), each = nSample)
batch.facLoc <- runif(nSample, min = 0, max = 0.5)
batch.facLoc <- rep(batch.facLoc, times =2)

################################################################################
############################## Useful Functions ################################
################################################################################
generate_theta <- function(nSample, nTimepoints = 2, nCellType, mean, sd, num_Cov = 1){
  Simplex <- nCellType - 1
  sample_ids <- rep(1:nSample, each = nTimepoints) 
  Timepoint <- rep(c(0, 1), nSample)
  Response <- rep(c(0, 1), each = (nSample * nTimepoints) / 2)
  # interaction Term
  Timepoint_Response <- Timepoint * Response
  
  # Create design matrix for each sample/time combination
  X <- data.frame(Timepoint = Timepoint,
                  Response = Response,
                  Timepoint_Response = Timepoint_Response)
  rownames(X) <- paste(paste0("Sample", sample_ids), ifelse(Response == 0, "nonResponse", "Response"),
                       ifelse(Timepoint == 0, "t1", "t2"), sep = "_")
  
  
  # Generate Effect Size Matrix
  # gamma <- data.frame(K1 = c(0, ResponseEffect, TimeEffect),
  #                     K2 = c(0, ResponseEffect, -TimeEffect))
  # gamma <- cbind(gamma, matrix(0, nrow = 3, ncol = (Simplex - 2)))
  gamma <- vector(mode = "list")
  for (i in 1:Simplex) {
    gamma[[i]] <- diag(mvrnorm(n = num_Cov, mu = mean[,i], Sigma = sd[[i]]))
  }
  gamma <- do.call(cbind, gamma)
  # gamma <- rbind(1, gamma)
  colnames(gamma) <- c(paste0("K", 1:ncol(gamma)))
  rownames(gamma) <- c("Timepoint", "Response", "Timepoint_Response")
  
  mu <- t(t(as.matrix(gamma))  %*% t(as.matrix(X)))
  
  # simulate patient random effect
  psi <- mvrnorm(n = nSample, mu = rep(0, Simplex),  Sigma= diag(rep(1, Simplex)))
  
  # adding random effect eta = mu + psi
  for (i in 1:nSample) {
    idx_t1 <- 2 * (i - 1) + 1
    idx_t2 <- 2 * i
    
    mu[idx_t1, ] <- mu[idx_t1, ] + psi[i, ]
    mu[idx_t2, ] <- mu[idx_t2, ] + psi[i, ]
  }
  
  # transfer the proportion into simplex.
  # adding 0 to every sample/time for identifiability in the model
  eta <- cbind(mu, 0)
  colnames(eta) <- paste0("K", 1:ncol(eta))
  
  map_to_simplx <- function(x) {
    exp(x - log(sum(exp(x))))
  }
  
  theta <- t(apply(eta, 1, map_to_simplx))
  theta <- list(t1 = theta[grep("t1", rownames(theta)), ], t2 = theta[grep("t2", rownames(theta)), ])
  # return(list(theta = theta, gamma = gamma, psi = psi))
  return(list(theta = theta, gamma = gamma))
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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse/"
if(batch.rmEffect){save_batch <- "noBatch"} else{save_batch <- "Batch"}
if(is.null(cancerCellGroup)){save_cancer <- "StromalCell"} else{save_cancer <- "CancerCell"}

r.file <- paste0("doc/rev_splatter/",list.files("doc/rev_splatter/"))
sapply(r.file, source)

param_dat <- expand.grid(numCellType, gamma_sd)
colnames(param_dat) <- c("numCellType", "gamma_sd")

# #############################################################
for (i in 1:nrow(param_dat)){
  nCellType <- param_dat[i,1]
  interaction_effect <- param_dat[i,2]
  simplex <- nCellType- 1 
  mean = matrix(rep(0, simplex*A), nrow = A)
  sd <- replicate(simplex, diag(c(rep(0, A - 1), interaction_effect), A), simplify = FALSE)
  if(interaction_effect == 0){
    type = "NullModel"
  }else{
    type = paste0("HighVar", interaction_effect)
  }
  
  
  true_param <- generate_theta(nSample = nSample, nTimepoints = 2, nCellType = nCellType, mean = mean, sd = sd, num_Cov = A)
  sims <- delta_sim(true_param = true_param, seed = seed, nSample = nSample,
                    nGenes = nGenes, nCellType = nCellType,
                    de.prob = de.prob, de.facLoc = de.facLoc,
                    batchCells = batchCells, batch.facLoc = batch.facLoc,
                    cancerCellGroup = cancerCellGroup,
                    batch.rmEffect = batch.rmEffect)
  
  # save_response <- paste(round(ResponseEffect, 1), collapse = ".")
  # save_time <- paste(round(TimeEffect, 1), collapse = ".")
  dir_path = paste0(dir, "nSample", nSample,
                    "_nCellType", nCellType, "_", save_batch, "_", save_cancer, "/sims/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  saveRDS(sims, file = paste0(dir_path,
                              "sims_", seed, "_", type, ".rds"))
  # print(proportion_check(sims))
  # print((true_param$theta$t2-true_param$theta$t1)/true_param$theta$t1)
  rm(sims)
  cat("Cell Type is", param_dat[i,1], "and the sd is", param_dat[i,2], "\n")
  # sims <- sims[,sims$Group == "Group2"]
  # sims <- logNormCounts(sims)
  # sims <- runPCA(sims)
  # plotPCA(sims, colour_by = "Group", shape_by = "Sample")
  # plotPCA(sims, colour_by = "Sample", shape_by = "Batch")
  # plotPCA(sims, colour_by = "Group")
  
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

