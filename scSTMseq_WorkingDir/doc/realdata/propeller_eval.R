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
nsims <- 100
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
collaposed.theta <- lapply(SimEta, function(x){
  x <- cbind(x,0)
  theta.old <- exp(x - log(rowSums(exp(x))))
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

collaposed.theta <- readRDS( "res/PD1/collaposed_theta.rds")
############################################################
##################### propeller ############################
############################################################
collaposed.theta <- readRDS( "res/PD1/collaposed_theta.rds")
TransformReplicates <- function(TopicProp, transform="logit"){
  props <- TopicProp %>%
    dplyr::select(starts_with("V"))
  if(transform=="asin"){
    message("Performing arcsin square root transformation of proportions")
    prop.trans <- asin(sqrt(props))
  }else if(transform=="logit"){
    message("Performing logit transformation of proportions")
    prop.trans <- log(props/(1-props))
  }
  colnames(prop.trans) <- paste0("Cluster", 1:ncol(prop.trans))
  prop.trans$Sample <- TopicProp$patient_id
  prop.trans$Time <- TopicProp$timepoint
  prop.trans$Response <- TopicProp$expansion
  # prop.trans <- t(prop.trans)
  return(prop.trans)
}

propellerReplicatesBlock <- function(prop.trans, trend=FALSE, robust=TRUE){
  sample <- prop.trans$Sample
  group <- factor(prop.trans$Time, levels = c("Pre", "On"))
  response <- factor(prop.trans$Response, levels = c("NE", "E"))
  replicates <- factor(prop.trans$Replicates)
  
  x.logit <- prop.trans %>% dplyr::select(starts_with("Cluster"))
  x.logit <- t(x.logit)
  cov <- factor(paste(response,group,sep="."))
  design<-model.matrix(~ 0 + cov)
  colnames(design) <- c("E_t2", "E_t1", "NE_t2", "NE_t1")
  
  corfit<-duplicateCorrelation(x.logit, design, block=replicates)
  fit<-lmFit(x.logit, design, block=replicates, correlation=corfit$consensus)
  
  contrasts <- makeContrasts(Response_Change = E_t2 - E_t1,
                             nonResponse_Change = NE_t2 - NE_t1,
                             Interaction = (E_t2 - E_t1) - (NE_t2 - NE_t1),
                             levels = design)
  fit.cont <- contrasts.fit(fit, contrasts = contrasts) 
  fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)
  
  out <- data.frame(coefficients = fit.cont$coefficients,
                    t = fit.cont$t,
                    df.total = fit.cont$df.total,
                    p.value = fit.cont$p.value)
  rownames(out) <- rownames(fit.cont$coefficients)
  out$fdr.Response_Change <- p.adjust(out$p.value.Response_Change, method="BH")
  out$fdr.nonResponse_Change <- p.adjust(out$p.value.nonResponse_Change, method="BH")
  out$fdr.Interaction <- p.adjust(out$p.value.Interaction, method="BH")
  return(out)
}


propRep.logit <- lapply(collaposed.theta, function(x){
  TransformReplicates(x, "logit")
})
names(propRep.logit) <- paste0("Replicates", 1:length(propRep.logit))
propRep.logit <- lapply(names(propRep.logit), function(name){
  x <- propRep.logit[[name]]
  x$Replicates <- name
  return(x)
})
propRep.logit <- do.call(rbind, propRep.logit)
res.logit <- propellerReplicatesBlock(propRep.logit)

saveRDS(res.logit, "res/PD1/propeller_BlockRep_logit.rds")

propRep.asin <- lapply(collaposed.theta, function(x){
  TransformReplicates(x, "asin")
})
names(propRep.asin) <- paste0("Replicates", 1:length(propRep.asin))
propRep.asin <- lapply(names(propRep.asin), function(name){
  x <- propRep.asin[[name]]
  x$Replicates <- name
  return(x)
})
propRep.asin <- do.call(rbind, propRep.asin)
res.asin <- propellerReplicatesBlock(propRep.asin)

saveRDS(res.logit, "res/PD1/propeller_BlockRep_asin.rds")
