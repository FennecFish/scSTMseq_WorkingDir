setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)
library(MASS)
library(tibble)
library(tidyverse)
library(MANOVA.RM)
library(limma)
source("R/useful_functions.R")
scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_BatchNoInteraction_Prevalence_TimeResponseInteraction.rds")
scSTMobj <- select_top_scSTM(scSTMobj)
nsims <- 1000
#########################################
## Step 1: simulate eta 
#########################################
source("R/etaPosterior.R")
model <- scSTMobj
SimEta <- etaPosterior(model = model, nsims = nsims)
metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, patient_id, timepoint, expansion)
#########################################
## Step 2: get topic proportion 
## then collapse by time and sample
## renormalize theta to proportion
#########################################
collapsed.theta <- lapply(SimEta, function(x){
  x <- cbind(x,0)
  max_x <- apply(x, 1, max) 
  stabilized_x <- sweep(x, 1, max_x, "-")
  
  theta.old <- exp(stabilized_x - log(rowSums(exp(stabilized_x))))
  rownames(theta.old) <- model$DocName
  theta.old <- theta.old %>%
    as.data.frame() %>%
    rownames_to_column("Cell") %>%
    left_join(metadata, by = "Cell")
  
  theta.collapsed <- theta.old %>%
    group_by(patient_id, timepoint, expansion) %>%
    summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
    ungroup()
  
  theta.new <- theta.collapsed %>%
    rowwise() %>%
    mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
    ungroup()
  
  return(theta.new)
})
saveRDS(collapsed.theta, file = "res/PD1/collaposed_theta.rds")

#########################################
## Step 3b: MANOVA + Repeated Measure
#########################################
# fit a manova model
collapsed.theta <- readRDS("res/PD1/collaposed_theta.rds")
manova.fit <- lapply(collapsed.theta, function(x){
  Y <- x %>% 
    dplyr::select(starts_with("V"))
  Y <- compositions::ilr(Y)
  x <- x %>%
    dplyr::select(-starts_with("V")) %>%
    dplyr::mutate(Time = factor(timepoint, levels = c("Pre", "On")),
                  Sample = as.factor(patient_id),
                  Response = factor(expansion, levels = c("NE", "E"))) %>%
    dplyr::select(Time, Sample, Response)
 #  x$Time <- as.numeric(x$Time)-1
  # x$Response <- as.numeric(x$Response)-1
  x <- cbind(x, Y)
  response_vars <- grep("^V", colnames(x), value = TRUE)
  warning_message <- "No Warnings"
  MANOVA.RM.res <- tryCatch({
    withCallingHandlers({
      multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Response*Time")),
             data = x, subject = "Sample", within = "Time", iter = 1000)},
      warning = function(w) {
        warning_message <<- conditionMessage(w)  # Capture the warning message
        invokeRestart("muffleWarning")  # Mute the warning
      })
  })
  MANOVA.RM.res$Warnings <- warning_message
  
  # MANOVA.res <- tryCatch({
  #   withCallingHandlers({
  #     MANOVA.wide(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Response*Time")),
  #                 data = x, iter = 1000)},
  #     warning = function(w) {
  #       warning_message <<- conditionMessage(w)  # Capture the warning message
  #       invokeRestart("muffleWarning")  # Mute the warning
  #     })
  # })
  return(MANOVA.RM.res)
})

saveRDS(manova.fit, file = "res/PD1/Manova_Res.rds")

