setwd("/proj/milovelab/wu/scLDAseq")
source("R/etaPosterior.R")

PosteriorReplicates <- function(model, nsims = 100){
  
  #########################################
  ## Step 1: simulate eta 
  #########################################
  SimEta <- etaPosterior(model, nsims = nsims)
  metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, Sample, Time, Response)
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
      group_by(Sample, Time) %>%
      summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
      ungroup()
    
    theta.new <- theta.collapsed %>%
      rowwise() %>%
      mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
      ungroup()
    
    return(theta.new)
  })
  
  return(collaposed.theta)
}

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
  prop.trans$Sample <- TopicProp$Sample
  prop.trans$Time <- TopicProp$Time
  # prop.trans <- t(prop.trans)
  return(prop.trans)
}

propellerReplicatesBlock <- function(prop.trans, trend=FALSE, robust=TRUE){
  sample <- prop.trans$Sample
  group <- factor(prop.trans$Time)
  replicates <- factor(prop.trans$Replicates)
  
  x.logit <- prop.trans %>% dplyr::select(starts_with("Cluster"))
  x.logit <- t(x.logit)
  design<-model.matrix(~0+factor(group))
  colnames(design) <- c("Time1", "Time2")
  
  corfit<-duplicateCorrelation(x.logit, design, block=replicates)
  fit<-lmFit(x.logit, design, block=replicates, correlation=corfit$consensus)
  
  contrasts <- makeContrasts(Time = Time2 - Time1, levels = design)
  fit.cont <- contrasts.fit(fit, contrasts = contrasts) 
  fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)
  
  out <- data.frame(coefficients = fit.cont$coefficients[,1],
                    t = fit.cont$t[,1],
                    df.total = fit.cont$df.total,
                    p.value = fit.cont$p.value[,1])
  rownames(out) <- rownames(fit.cont$coefficients)
  out$fdr <- p.adjust(out$p.value, method="BH")
  return(out)
}

propellerReplicates <- function(prop.trans, trend=FALSE, robust=TRUE){
  sample <- prop.trans$Sample
  group <- factor(prop.trans$Time)
  design<-model.matrix(~0+factor(group))
  colnames(design) <- c("Time1", "Time2")
  # corfit<-duplicateCorrelation(prop.trans,design,block=sample)
  props <- prop.trans %>% dplyr::select(starts_with("Cluster"))
  props <- t(props)
  fit<-lmFit(props,design)
  contrasts <- makeContrasts(Time = Time2 - Time1, levels = design)
  fit.cont <- contrasts.fit(fit, contrasts = contrasts) 
  fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)
  out <- data.frame(coefficients = fit.cont$coefficients[,1],
                    t = fit.cont$t[,1],
                    df.total = fit.cont$df.total,
                    p.value = fit.cont$p.value[,1])
  rownames(out) <- rownames(fit.cont$coefficients)
  out$fdr <- p.adjust(out$p.value, method="BH")
  return(out)
}
