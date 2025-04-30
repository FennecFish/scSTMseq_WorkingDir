# This script is to estimate measurement error in Gamma

GammaError <- function(model, formula, null_formula = NULL, nsims = 100){
  # simulate eta
  SimEta <- etaPosterior(model, nsims = nsims)
  metadata <- colData(model$settings$sce)
  
  metadata <- colData(model$settings$sce)
 
  ##
  #Step 1: Extract the formula and do some error checking
  ##
  if(!inherits(formula,"formula")) stop("formula must be a formula object.")
  if(!is.null(metadata) & !is.data.frame(metadata)) metadata <- as.data.frame(metadata)
  
  # if (!is.null(ref.vec)) {
  #   # Extract the covariate terms from the formula
  #   covariate_terms <- strsplit(as.character(formula)[3], "\\s*\\*\\s*")[[1]]
  #   if (length(covariate_terms) != length(ref.vec)) stop("The size of the reference group should be the same as the number of covariates.")
  #   
  #   # Apply reference and alternative to each covariate and reference pair
  #   for (i in seq_along(covariate_terms)) {
  #     covariate <- covariate_terms[i]
  #     if(!covariate %in% names(ref.vec)) stop("The name of the reference group vector should be the covariate specified.")
  #     ref_value <- ref.vec[match(covariate, names(ref.vec))]
  #     if (!ref_value %in% metadata[[covariate]]) stop("The reference value '", ref_value, "' should be an element in the covariate '", covariate, "'")
  #     # Update the factor levels with the reference value as the first level
  #     metadata[[covariate]] <- factor(metadata[[covariate]], levels = c(ref_value, setdiff(unique(metadata[[covariate]]), ref_value)))
  #   }
  # }
  
  # check for random effect
  terms <- unlist(strsplit(deparse(formula), "\\+")) 
  random_effects <- gsub("\\s+", "", terms[grepl("\\|", terms)])
  random_effects <- gsub(".*\\|(.*)\\)", "\\1", random_effects)
  
  if(length(random_effects) > 0){
    fixed_effects <- terms[!grepl("\\|", terms)]
    fixed_effects <- if (length(fixed_effects) > 0) {
      as.formula(paste(paste(fixed_effects, collapse = " + ")))

    } else {
      as.formula("~ 1")  # If there are no fixed effects, return an intercept-only formula
    }
    termobj <- terms(fixed_effects, data=metadata)
    xmat <- try(Matrix::sparse.model.matrix(termobj,data=metadata),silent=TRUE)
    xmat <- as.matrix(xmat)
    # if random effect exists, fit a mixed linear model
    mixed.covar <- xmat %>%
      as.data.frame() %>%
      mutate(Sample = metadata$Sample[match(rownames(xmat), metadata$Cell)])
    X <- mixed.covar[,-1]
    gamma <- vector(mode = "list")
    lm.model <- vector(mode = "list")
    null.model <- vector(mode = "list")
    model.diff <- vector(mode = "list")
    random_effect <- paste0("(1|", random_effects, ")")
    lmer.formula <- paste(c(colnames(xmat)[-1], random_effect), collapse = "+")
    lmer.formula <- as.formula(paste0("Y_all[, k]~", lmer.formula))
    
    if(!is.null(null_formula)){
      null_term <- unlist(strsplit(deparse(null_formula), "\\+")) 
      null_fixed_effects <- null_term[!grepl("\\|", null_term)]
      if(length(null_fixed_effects) >0){
        null_fixed_effects <- as.formula(paste(null_fixed_effects, collapse = " + "))
        termobj.null <- terms(null_fixed_effects, data=metadata)
        xmat.null <- try(Matrix::sparse.model.matrix(termobj.null,data=metadata),silent=TRUE)
        xmat.null <- as.matrix(xmat.null)
        # null.mixed.covar <- xmat.null %>%
        #   as.data.frame() %>%
        #   mutate(Sample = metadata$Sample[match(rownames(xmat.null), metadata$Cell)])
        # X.null <- null.mixed.covar[,-1]
        null.formula <- paste(c(colnames(xmat.null)[-1], random_effect), collapse = "+")
        null.formula <- as.formula(paste0("Y_all[, k]~", null.formula))
      }else{
        # X.null <- mixed.covar %>% dplyr::select(random_effects)
        # X.null <- X
        null.formula <- as.formula(paste0("Y_all[, k]~", random_effect))
      }
    }

    pb <- txtProgressBar(min = 0, max = length(SimEta), style = 3)
    for(i in 1:length(SimEta)){
      Y_all <- SimEta[[i]] 
      res <- data.frame()
      lm.model.k <- vector(mode = "list")
      null.model.k <- vector(mode = "list")
      model.diff.k <- vector(mode = "list")
      for (k in seq_len(ncol(Y_all))) {
        lm.model.temp <- lme4::lmer(lmer.formula, data = X)
        temp.res <- summary(lm.model.temp)$coefficients[,"Estimate"]
        res <- rbind(res, temp.res)
        lm.model.k[[k]] <- lm.model.temp
        null.model.k[[k]] <- lme4::lmer(null.formula, data = X)
        model.diff.k[[k]] <- anova(lm.model.temp, null.model.k[[k]])
      }
      colnames(res) <- names(temp.res)
      rownames(res) <- paste0("topic_", 1:ncol(Y_all))
      gamma[[i]] <-  res
      lm.model[[i]] <- lm.model.k
      null.model[[i]] <- null.model.k
      model.diff[[i]] <- model.diff.k
      setTxtProgressBar(pb, i)
    }
}
  toreturn <- list(SimEta = SimEta, EstGamma = gamma, formula = formula, full.model = lm.model, null.model = null.model,
                   model.diff = model.diff)
  return(toreturn)
}
