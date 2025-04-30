# Goal: This script is to generate cell type proportions before- and after- treatment
#       given different response group, for multiple patients
#       Response group will have fluctuation in cell type proportion after treatment
#       non-response group will NOT have fluctuation in cell type proportion after treatment
#       Both pre and post treatment proportion is generated from a Logistic Normal, 
#       followed by softmax transformation
#       We then visualize this difference between response group
# Author: Euphy Wu
# Date: Oct 30th, 2025
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(tidyverse)
library(MASS)
library(stringr)
################################################################################
############################## Useful Functions ################################
################################################################################
# This function is used to generate cell type proportion for before- and after- treatment
# nSample: number of patients
# nTimepoints: number of timepoints, defult = 2
# nCellType: Number of Cell types
# mean: mean of effect size
# sd: variance of effect size
# num_Cov:number of covariate of interest
generate_theta <- function(nSample, nTimepoints = 2, nCellType, mean, sd, num_Cov = 1){
    # To simulate proportion data, first we work in nCellType - 1 simplex
    Simplex <- nCellType - 1
    # define design matrix
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
    
    # Generate Effect Size Matrix, by sampling from Multivariate normal distribution
    # given pre-defined mean and covariance matrix
    gamma <- vector(mode = "list")
    for (i in 1:Simplex) {
        gamma[[i]] <- diag(mvrnorm(n = num_Cov, mu = mean[,i], Sigma = sd[[i]]))
    }
    gamma <- do.call(cbind, gamma)
    colnames(gamma) <- c(paste0("K", 1:ncol(gamma)))
    rownames(gamma) <- c("Timepoint", "Response", "Timepoint:Response")
    
    # Generate mu = X*effect size
    mu <- t(t(as.matrix(gamma))  %*% t(as.matrix(X)))
    
    # simulate patient random effect psi from a standard normal
    psi <- mvrnorm(n = nSample, mu = rep(0, Simplex),  Sigma= diag(rep(1, Simplex)))
    
    # adding random effect to the pre-calculated mean: eta = mu + psi
    for (i in 1:nSample) {
        idx_t1 <- 2 * (i - 1) + 1
        idx_t2 <- 2 * i
        
        mu[idx_t1, ] <- mu[idx_t1, ] + psi[i, ]
        mu[idx_t2, ] <- mu[idx_t2, ] + psi[i, ]
    }
    
    # softmax transformation of the simplex to proportion
    # adding 0 to every sample/time for identifiability in the model
    eta <- cbind(mu, 0)
    colnames(eta) <- paste0("K", 1:ncol(eta))
    
    map_to_simplx <- function(x) {
        exp(x - log(sum(exp(x))))
    }
    
    theta <- t(apply(eta, 1, map_to_simplx))
    theta <- list(t1 = theta[grep("t1", rownames(theta)), ], t2 = theta[grep("t2", rownames(theta)), ])
    return(list(theta = theta, gamma = gamma))
}

###############################################################################
########################### Simulation ########################################
###############################################################################

################## Define Parameters ###################
seed = 123
numCellType <- 5
# the variance of logistic normal dist.
# we want to test the effect given different effect size
gamma_sd <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# number of patients
nSample <- 20
# since we are interested in before- and after- treatment
# timepoint is 2
nTimepoints <- 2

# number of covariates in the logistic normal regression
# Response, Time and interaction between response:time
A = 3 

param_dat <- expand.grid(numCellType, gamma_sd)
colnames(param_dat) <- c("numCellType", "gamma_sd")

################### Data Generation ##########################
true_param <- vector(mode = "list")
# generate corresponding dataset given different covariance/ number of cell type combo
for (i in 1:nrow(param_dat)){
  nCellType <- param_dat[i,1]
  interaction_effect <- param_dat[i,2]
  simplex <- nCellType- 1 
  mean = matrix(rep(0, simplex*A), nrow = A)
  # since nonresponse group does not have cell type proportion change
  # the effect size for both time and interaction should be 0
  sd <- replicate(simplex, diag(c(rep(0, A - 1), interaction_effect), A), simplify = FALSE)
  true_param[[i]] <- generate_theta(nSample = nSample, nTimepoints = 2, nCellType = nCellType, mean = mean, sd = sd, num_Cov = A)
}

# For each parameter combination, we want to calculate the relative change in cell type proportion 
# averaged across samples within each group
proportion_change <- lapply(true_param, function(x){
    # calculate relative change
  RC <- (x$theta$t2-x$theta$t1)/x$theta$t1 %>% 
    as.data.frame()
  # extract response group
  RC$Response =str_extract(rownames(RC), "(?<=_)[^_]+(?=_)")
  RC <- RC %>%
    group_by(Response) %>%
    summarise(across(starts_with("K"), ~mean(.x, na.rm = TRUE)))
})
names(proportion_change) <- gamma_sd

###############################################################################
############################### Plot ##########################################
###############################################################################
# row bind the list and then pivot to long format
proportion_change <- bind_rows(proportion_change, .id = "gamma") %>%
  mutate(gamma = as.numeric(gamma)) %>%
  pivot_longer(cols = starts_with("K"), names_to = "K", values_to = "value")
proportion_change$gamma <- paste0("Effect Size = ", proportion_change$gamma)

ggplot(proportion_change, aes(x = K, y = value, color = Response, group = Response)) +
    geom_line() +
    ggtitle("Relative Change in Cell Type Proportion Across Different Response Groups and Effect Size") +
    facet_wrap(~ gamma) +  
    labs(x = "Cell Type", y = "Relative Change", color = "Response Group") +
    theme_bw() +
    theme(strip.text = element_text(size = 12), legend.position = "bottom")
 