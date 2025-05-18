# This script is to apply manova with repeated measurement
library(MASS)
library(MANOVA.RM)
# library(CompDTUReg)
library(stats)
library(dplyr)
library(SingleCellExperiment)
library(tibble)
# library(MCMCglmm)

# simualte data and collapose by Sample and temporal information
PosteriorPropRep <- function(model, nsims = 100, Sample = "Sample", Time = "Time", Group = NULL){
    
  #########################################
  ## Step 1: simulate eta 
  #########################################
  SimEta <- etaPosterior(model, nsims = nsims)

  if(!is.null(Group)){
    required_cols <- c("Cell", Sample, Time, Group)
    metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(all_of(required_cols))
    #########################################
    ## Step 2: get topic proportion 
    ## then collapse by time and sample
    ## renormalize theta to proportion
    #########################################
    collaposed.theta <- lapply(SimEta, function(x){
      x <- cbind(x,0)
      # to avoid exp(x) leading to inf, we stablize the value by subtracting the max
      max_x <- apply(x, 1, max) 
      stabilized_x <- sweep(x, 1, max_x, "-")
      
      theta.old <- exp(stabilized_x - log(rowSums(exp(stabilized_x))))
      rownames(theta.old) <- model$DocName
      theta.old <- theta.old %>%
        as.data.frame() %>%
        rownames_to_column("Cell") %>%
        left_join(metadata, by = "Cell")
      
      theta.collapsed <- theta.old %>%
        group_by(across(all_of(c(Sample, Time, Group)))) %>%
        summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
        ungroup()
      
      theta.new <- theta.collapsed %>%
        rowwise() %>%
        mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
        ungroup()
      
      return(theta.new)
    })
  }else{
    required_cols <- c("Cell", Sample, Time)
    metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(all_of(required_cols))
    #########################################
    ## Step 2: get topic proportion 
    ## then collapse by time and sample
    ## renormalize theta to proportion
    #########################################
    collaposed.theta <- lapply(SimEta, function(x){
      x <- cbind(x,0)
      # to avoid exp(x) leading to inf, we stablize the value by subtracting the max
      max_x <- apply(x, 1, max) 
      stabilized_x <- sweep(x, 1, max_x, "-")
      
      theta.old <- exp(stabilized_x - log(rowSums(exp(stabilized_x))))
      rownames(theta.old) <- model$DocName
      theta.old <- theta.old %>%
        as.data.frame() %>%
        rownames_to_column("Cell") %>%
        left_join(metadata, by = "Cell")
      
      theta.collapsed <- theta.old %>%
        group_by(across(all_of(c(Sample, Time)))) %>%
        summarise(across(starts_with("V"), sum), .groups = 'drop') %>%
        ungroup()
      
      theta.new <- theta.collapsed %>%
        rowwise() %>%
        mutate(across(starts_with("V"), ~ . / sum(c_across(starts_with("V"))))) %>%
        ungroup()
      
      return(theta.new)
    })
  }

  return(collaposed.theta)
}

# fit manova with posterior replicates
ThetaManova <- function(PosteriorRep, Time = "Time", Group = NULL, use_mean = T){
  if(use_mean){
    combined <- do.call(rbind, PosteriorRep)
    PosteriorRep <- combined %>% 
      group_by(Sample, Time) %>%
      summarise(across(starts_with("V"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    PosteriorRep <- list(PosteriorRep)
  }
  #########################################
  ## Step 3b: MANOVA + Repeated Measure
  #########################################
  # fit a manova model
  manova.fit <- lapply(PosteriorRep, function(data){
    Y <- data %>% 
      dplyr::select(starts_with("V"))
    Y <- compositions::ilr(Y)
    data <- data %>%
      dplyr::select(-starts_with("V")) %>%
      dplyr::mutate(Time = as.factor(Time),
             Sample = as.factor(Sample))
    data <- cbind(data, Y)
    response_vars <- grep("^V", names(data), value = TRUE)
    warning_message <- "No Warnings"
    MANOVA.RM.res <- tryCatch({
        withCallingHandlers({
          if(is.null(Group)){covariate <- Time}
          else{covariate <- paste(c(Group, Time), collapse = "*")}
          multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~", covariate)),
                 data = data, subject = "Sample", within = "Time", iter = 1000)
          },
          warning = function(w) {
            warning_message <<- conditionMessage(w)  # Capture the warning message
            invokeRestart("muffleWarning")  # Mute the warning
          })
      })
    MANOVA.RM.res$Warnings <- warning_message
    
    # check if interaction term is significant
    # if not significant, run main effect only
    interaction_term <- grepl(":", rownames(MANOVA.RM.res$resampling))
    if(sum(interaction_term) > 0){
      interact.p <- MANOVA.RM.res$resampling[interaction_term, 2]}else{
        interact.p <- 0
      }
    if(!is.null(Group) & interact.p > 0.05){
      warning_message <- "No Warnings"
      MANOVA.RM.main <- tryCatch({
        withCallingHandlers({
          covariate <- paste(c(Group, Time), collapse = "+")
          multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~", covariate)),
                 data = data, subject = "Sample", within = "Time", iter = 1000)
        },
        warning = function(w) {
          warning_message <<- conditionMessage(w)  # Capture the warning message
          invokeRestart("muffleWarning")  # Mute the warning
        })
      })
      MANOVA.RM.main$Warnings <- warning_message
    }else{
      MANOVA.RM.main <- NA
    }
    # MANOVA.res <- tryCatch({
    #   withCallingHandlers({
    #     if(is.null(Group)){covariate <- Time }
    #     else{covariate <- paste(c(Group, Time), collapse = "*")}
    #     MANOVA.wide(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~", covariate)),
    #                 data = data, iter = 1000)},
    #     warning = function(w) {
    #       warning_message <<- conditionMessage(w)  # Capture the warning message
    #       invokeRestart("muffleWarning")  # Mute the warning
    #     })
    # })
    return(list(Interaction = MANOVA.RM.res, MainEffect = MANOVA.RM.main))
  })
  
  return(manova.fit)
}

MaxThetaManova <- function(model, nsims = 100, Time = "Time", Sample = "Sample"){
  
  SimEta <- etaPosterior(model, nsims = nsims)
  
  max_cluster <- lapply(SimEta, function(rep){
    max_indices <- apply(rep, 1, which.max)
    colnames(rep) <- paste0("topic_", 1:ncol(rep))
    rownames(rep) <- model$DocName
    res_cluster <- colnames(rep)[max_indices]
    names(res_cluster) <- rownames(rep)
    
    metadata <- colData(model$settings$sce) %>% as.data.frame()
    metadata$Cluster <- res_cluster[match(rownames(metadata), names(res_cluster))]
    
    max_proportions <- metadata %>%
      group_by(Sample, Time, Cluster) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(Sample, Time) %>%
      mutate(proportion = count / sum(count)) %>%
      ungroup() %>%
      dplyr::select(-count) %>%
      pivot_wider(names_from = Cluster, values_from = proportion, values_fill = 0)
    
    meta <- max_proportions %>% dplyr::select(!starts_with("topic"))
    cluster <- max_proportions %>% dplyr::select(starts_with("topic"))
    
    cluster <- compositions::ilr(cluster)
    colnames(cluster) <- paste0("topic", 1:ncol(cluster))
    data <- cbind(meta, cluster)
    ##### Fit Manova.RM #######
    response_vars <- grep("^topic", colnames(data), value = TRUE)
    warning_message <- "No Warnings"
    fit <- tryCatch({
      withCallingHandlers({
        covariate <- Time
        multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~", covariate)),
               data = data, subject = "Sample", within = "Time", iter = 1000)
      },
      warning = function(w) {
        warning_message <<- conditionMessage(w)  # Capture the warning message
        invokeRestart("muffleWarning")  # Mute the warning
      })
    })
    fit$Warnings <- warning_message
    return(fit)
  })
  
  return(max_cluster)
}
  # #########################################
  # ## Step 3c: multivariate mixed linear effect
  # #########################################
  # # fit a manova model
  # mcmcglmm.fit <- lapply(collaposed.theta, function(x){
  #   Y <- x %>% 
  #     dplyr::select(starts_with("V"))
  #   Y <- compositions::ilr(Y)
  #   x <- x %>%
  #     dplyr::select(-starts_with("V")) %>%
  #     dplyr::mutate(Time = as.factor(Time),
  #                   Sample = as.factor(Sample))
  #   x <- cbind(x, Y)
  #   response_vars <- grep("^V", names(x), value = TRUE)
  #   fixed_formula <- as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time "))
  #   ilr_dim <- ncol(Y)
  #   priors = list(R = list(V = diag(ilr_dim), nu = 0.002),
  #                 G = list(G1 = list(V = diag(ilr_dim), nu = 2,
  #                                    alpha.mu = rep(0, ilr_dim),
  #                                    alpha.V = diag(1, ilr_dim, ilr_dim))))
  #   mcmc.fit <- tryCatch(
  #     MCMCglmm(fixed_formula, random = ~us(trait):Sample, rcov = ~us(trait):units, data = x,
  #       family = rep("gaussian", ilr_dim), verbose = FALSE, prior = priors),
  #     error = function(e) NA)
  # })
  # 


# This function is to calculate pValue directly copied from compDTU
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, lm_model_fit = NULL, q, nsamp, df_residual = NA){
  
  
  # #########################################
  # ## Step 3a: Directly using CompDTU
  # #########################################
  # # calculate irl for each replicate
  # irl.InfRep <- lapply(1:length(collaposed.theta), function(i){
  #   y <- collaposed.theta[[i]] %>% dplyr::mutate(InfRep = paste0("infRep", i))
  #   Y <- y %>% 
  #     dplyr::select(starts_with("V"))
  #   Y <- compositions::ilr(Y)
  #   res <- y %>% 
  #     dplyr::select(-starts_with("V"))
  #   res <- cbind(res, Y)
  #   return(res)
  # })
  # irl.InfRep <- do.call(rbind, irl.InfRep)
  # 
  # # calculate within sample covariance for each sample
  # WithinSampleCov <- vector(mode = "list")
  # unique.sample <- unique(irl.InfRep$Sample)
  # Z <- stats::model.matrix(~unique.sample)
  # for (sample in unique.sample){
  #   irl.sample <- irl.InfRep %>% 
  #     dplyr::filter(Sample == sample) %>%
  #     dplyr::select(starts_with("V"))
  #   ilrCov <- stats::cov(irl.sample)
  #   WithinSampleCov[[sample]] <- ilrCov
  # }
  # # calculate mean across samples
  # WithinSampleCov <- Reduce(`+`, WithinSampleCov) / length(WithinSampleCov)
  # TimeEffect <- irl.InfRep$Time
  # ncond <- length(unique(TimeEffect))
  # nsamp <- length(unique.sample)
  # # the following code are adapted from CompDTU
  # 
  # # alternative model
  # XAlt <- stats::model.matrix(~TimeEffect)
  # YInfRep <- irl.InfRep %>% dplyr::select(starts_with("V")) %>% as.matrix()
  # ns <- nrow(YInfRep)
  # XAltT <- t(XAlt)
  # bhatalt <- solve(crossprod(XAlt)) %*% (XAltT %*% YInfRep)
  # pie1 <- YInfRep - (XAlt%*%bhatalt)
  # SigmaTildeAltNewModeling <- (crossprod(pie1))/ns
  # 
  # # null model
  # XNull <- stats::model.matrix(~1, data = irl.InfRep)
  # XNullT <- t(XNull)
  # bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% YInfRep)
  # pie2 <- YInfRep - (XNull%*%bhatnull)
  # SigmaTildeNullNewModeling <- (crossprod(pie2))/ns
  # 
  # # Between_Sample_covariance = covariance - Within Sample Covariance
  # UpdatedCovAlt <- SigmaTildeAltNewModeling - WithinSampleCov
  # UpdatedCovNull <- SigmaTildeNullNewModeling - WithinSampleCov
  # 
  # #If there is a negative variance term, the pvalue for CompDTUme will be undefined
  # statement1 <- sum(diag(UpdatedCovAlt)<=0) !=0
  # statement2 <- sum(diag(UpdatedCovNull)<=0) !=0
  # ret.mixed <- data.frame(pval_CompDTUme = NA, FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)
  # ret <- data.frame(pval_CompDTUme = NA, FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)
  # 
  # if((statement1==TRUE | statement2==TRUE)){
  #   ret.mixed$pval_CompDTUme <- NA
  #   ret$pval_CompDTUme <- NA
  # }else{
  #   qq <- ncond - 1
  #   CompDTUmeRes.mixed <- tryCatch(calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
  #                                           lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
  #                            error = function(x){})
  #   
  #   if(is.null(CompDTUmeRes.mixed)==FALSE){
  #     ret.mixed$pval_CompDTUme <- CompDTUmeRes.mixed$pval_pillai
  #     ret.mixed$FStat <- CompDTUmeRes.mixed$FStat
  #     ret.mixed$NumDF <- CompDTUmeRes.mixed$NumDF
  #     ret.mixed$DenomDF <- CompDTUmeRes.mixed$DenomDF
  #   }
  #   CompDTUmeRes <- tryCatch(calcPillaiPval(SigmaTildeNull = SigmaTildeNullNewModeling, SigmaTildeAlt = SigmaTildeAltNewModeling,
  #                                                 lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
  #                                  error = function(x){})
  #   if(is.null(CompDTUmeRes)==FALSE){
  #     ret$pval_CompDTUme <- CompDTUmeRes$pval_pillai
  #     ret$FStat <- CompDTUmeRes$FStat
  #     ret$NumDF <- CompDTUmeRes$NumDF
  #     ret$DenomDF <- CompDTUmeRes$DenomDF
  #   }
  # }
  
  if(is.null(lm_model_fit) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in lm_model_fit to extract df.residual from")
  }
  
  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)
  # Etilde <- SigmaTildeAlt
  # Htilde <- SigmaTildeNull - SigmaTildeAlt
  
  vals <- diag((Htilde %*%solve(Etilde + Htilde)))
  pill_stat <- sum(vals)
  
  
  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas
  #This proved to be the easiest way to calculate the statistic and are confirmed to match R's anova.mlm function
  #v <- nsamp - ncol(stats::model.matrix(lm_model_fit))
  
  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- lm_model_fit$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit
  
  
  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)
  
  s <- min(p,q)
  
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)
  
  #Formulas come from the SAS help file Multivariate ANOVA Testing and are confirmed to match R's anova.mlm function
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))
  
  if(fstat_pillai < 0){
    return(list(FStat = NA, NumDF = NA, DenomDF = NA, pval_pillai = NA))
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-stats::pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  return(list(FStat = fstat_pillai, NumDF = numdf_pillai, DenomDF = denomdf_pillai, pval_pillai = pval_pillai))
}
