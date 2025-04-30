thetaPosteriorSample <- function(model, nsims=100) {
    lambda <- model$eta
    nu <- model$nu
    out <- vector(mode="list",length=nrow(lambda)) 
    for (i in 1:length(out)) {
        sigma <- nu[[i]]
        choleskydecomp <- chol(sigma)
        mat <- rmvnorm(nsims, lambda[i,],nu[[i]],choleskydecomp)
        mat <- cbind(mat, 0)
        out[[i]] <- exp(mat - row.lse(mat))
    }
    return(out)
}

mixed.lm <- function(stmobj, xmat, K, formula) {
  thetasims <- thetaPosteriorSample(stmobj, nsims=1) # draw from posterior distribution
  thetasims <- do.call(rbind, thetasims)
  # compute 
  thetaLogit <- log(thetasims/(1-thetasims))
  output <- vector(mode="list", length=length(K))
  # limma test
  
  for(k in K) {
    y <- thetaLogit[,k]
    # lm.mod = lmer(y ~ time + (1|sample), REML = F, data = xmat)
    
    form <- as.formula(paste0("y ~ ", as.character(formula)[-1], " + (1|sample)"))
    lm.mod = lmer(form, REML = F, data = xmat)
    # est <- summary(lm.mod)$coefficients[,1]
    lm_0 <-lmer(y ~ (1|sample), REML = F, data = xmat)
    anova(lm.mod,lm_0)
    PBmodcomp(lm.mod,lm_0,nsim=200)
    output[[which(k==K)]] <- est
  }
  names(output) <- paste0("topic",K)
  output <- do.call(rbind,output)
  return(output)
}

estimateMixedEffect <- function(formula, 
                           stmobj, metadata=NULL, 
                           # slope = FALSE,
                           # sampleNames = NULL, sampleIDs = NULL, # if id = null, then combine all samples
                           uncertainty=c("Global", "None"), 
                           nsims=25, numCores = TRUE) {
    origcall <- match.call()
    thetatype <- match.arg(uncertainty)
    
    if(thetatype=="None") nsims <- 1 #override nsims for no uncertainty
    
    #Step 1: Extract the fixed effect formula and do some error checking
    ##
    if(!inherits(formula,"formula")) stop("formula must be a formula object.")
    # if(!is.null(null.formula) & !inherits(null.formula,"formula")) stop("formula must be a formula object.")
    if(!is.null(metadata) & !is.data.frame(metadata)) metadata <- as.data.frame(metadata)
    if(is.null(metadata)) {
      metadata = as.data.frame(stmobj$settings$covariates$X)
      colnames(metadata) <- sub(".*\\$(.*)", "\\1", colnames(metadata))
    } 
    termobj <- terms(formula, data=metadata)
    if(attr(termobj, "response")==1){
        #if a response is specified we have to parse it and remove it.
        # as.character of a formula turns
        # dv ~ iv
        # into:  c("~", "dv", "iv")
        response <- as.character(formula)[2] #second object is the response in this cases
        K <- eval(parse(text=response))
        # browser()
        if(!(posint(K) && max(K)<=stmobj$settings$dim$K)) stop("Topics specified as response in formula must be a set of positive integers equal to or less than the number of topics in the model.")   
        #now we reconstruct the formula removing the response
        formula <- formula(paste(as.character(formula)[c(1,3)], collapse = " "))
 
        #the above used to be the below code but the use got deprecated.
        #as.formula(as.character(formula)[c(1,3)])
        termobj <- terms(formula, data=metadata)
    } else {
        K <- 1:stmobj$settings$dim$K
    }
    
    # now we subset the metadata to include only sample specified
    # if(!is.null(sampleIDs)) {
    #     if(is.null(sampleNames)) stop("Please specify the colname name for sampleIDs in the input metadata")
    #     if(!sampleIDs %in% stmobj$sampleID) stop("Sample ids specified must be exactly the same as the ones used in the model")
    #     metadata = metadata[metadata[,sampleNames] == sampleIDs,]
    #     subDocName <- stmobj$DocName[stmobj$sampleID %in% sampleIDs]
    #     stmobj <- STMsubset(stmobj, subDocName) 
    # } 
    mf <- model.frame(termobj, data=metadata)
    xmat <- model.matrix(termobj,data=metadata)
    varlist <- all.vars(termobj)
    # if(!is.null(metadata)) {
    #     data <- metadata[, varlist, drop=FALSE]
    # } else {
    #     templist <- list()
    #     for(i in 1:length(varlist)) {
    #         templist[[i]] <- get(varlist[i])
    #     }
    #     data <- data.frame(templist)
    #     names(data) <- varlist
    #     rm(templist)
    # }
    metadata <- metadata[, varlist, drop=FALSE]
    # rm(data)
    
    ##
    #Step 2: Draw Response Variables from Posterior Distribution
    ##
    xmat <- as.data.frame(xmat)
    xmat$sample <- stmobj$sampleID
    # xmat <- xmat[,-1]
    
    thetasims <- thetaPosteriorSample(stmobj, nsims=1) # draw from posterior distribution
    thetasims <- do.call(rbind, thetasims)
    # limma test
    thetasims <- t(thetasims)
    rownames(thetasims) <- paste0("topic",K)
    
    Time <- factor(xmat[,2])
    levels(Time) <- c("pre", "on")
    tdesign <- model.matrix(~0+Time)
    colnames(tdesign) <- levels(Time)
    corfit <- duplicateCorrelation(thetasims,tdesign,block=stmobj$sampleID)
    fit <- lmFit(thetasims,design = tdesign,block=stmobj$sampleID,correlation=corfit$consensus)
    cm <- makeContrasts(Time1v2 = on-pre, levels=tdesign)
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    topTable(fit2, coef="Time1v2")
    
    fit1$coefficients
    # target <- data.frame(Subject =  rep(1:6, each =2), 
    #                      +                      Condition = rep(c("Diseased", "Normal"), each = 6),
    #                      +                      Tissue = rep(c("A", "B"), rep = 6))
    # colnames(thetasims) <- pas
    # compute 
    # thetaLogit <- log(thetasims/(1-thetasims))
    # output <- vector(mode="list", length=length(K))
    
    
    ##  
    #Step 3: Calculate Coefficients
    ##
    
    # first simulate theta 
    storage <- vector(mode="list", length=nsims)
    
    # setup parallel running if numCores is specified
    if(numCores|is.numeric(numCores)) {
      if(is.numeric(numCores)) {
        cl <- numCores 
        if(cl > detectCores()) stop(paste0("Number of cores is more than available. The max cores available is ", detectCores()))
      } else{
        numCores <- detectCores()
      }
      cl <- makeCluster(numCores)
      clusterEvalQ(cl, {
        library(lme4)
      })
      clusterExport(cl, c("stmobj","xmat","K", "mixed.lm", "thetaPosteriorSample", "rmvnorm",
                          "row.lse", "formula"))
      # start_time <- Sys.time()
      storage <- parLapply(cl, 1:nsims,function(i) mixed.lm(stmobj, xmat, K, formula))
      stopCluster(cl)
      # end_time <- Sys.time()
      # elapsed_time <- end_time - start_time
    } else{
      
      for(i in 1:nsims) {
        storage[[i]] <- mixed.lm(stmobj, xmat, K)
        if (i %% 10 == 0) {  # Print progress every 10 iterations
          cat(sprintf("...%d%% ", i /nsims * 100))
          flush.console()
        }
      }
      cat("\n")
    }
    
    # storage <- do.call(rbind, storage)
    # transform theta into logit of theta
    
    # thetaLogit <- log(thetasims/(1-thetasims))
    
    # 3b) perform linear mixed effect model

    # browser()
    ##
    #Step 4: Return Values
    ##

    # Process 'output' for 'plrt'
    # plrt <- lapply(K, function(k) as.data.frame(output[[k]]$LRT))
    # names(plrt) <- paste0("topic", K)
    
    # Process 'output' for 'param'
    # param <- lapply(K, function(k) list(est = output[[k]]$coef[, 1], # Assuming you want all rows, first column for estimates
    #                                     std = output[[k]]$coef[, 2], # Assuming you want all rows, second column for std
    #                                     vcov = output[[k]]$vcov))
    # names(param) <- paste0("topic", K)
    # 
    # # Process 'output' for 'fixEffect'
    # fixEffect <- lapply(K, function(k) output[[k]]$coef)
    # names(fixEffect) <- paste0("topic", K)
    # 
    toreturn <- list(param.est=storage, 
                     topics=K,
                     data=metadata,
                     modelframe=mf, varlist=varlist)
    class(toreturn) <- "estimateMixedEffect"
    return(toreturn)
}

summary.estMixedEffect <- function(object, topics = NULL){
  if(!class(object)=="estimateMixedEffect")stop("The input needs to be a estimateMixedEffect Object")
  if(is.null(topics)) topics <- object$topics
  if(any(!(topics %in% object$topics))) {
    stop("Some topics specified with the topics argument are not available in this estimateMixedEffect object.")
  }
  tables <- vector(mode="list", length=length(topics))
  for(i in 1:length(topics)) {
    topic <- topics[i]
    est.sim <- lapply(object$param.est, function(df) df[paste0("topic", topic), ])
    est.sim <- do.call(rbind, est.sim)
    # apply t-test to each covariate
    output <- apply(est.sim, 2, function(x) t.test(x, mu = 0))
    est <- do.call(rbind,lapply(output, function(df) df["estimate"]))
    stdr <- do.call(rbind,lapply(output, function(df) df["stderr"]))
    tvalue <- do.call(rbind,lapply(output, function(df) df["statistic"]))
    p.value <- do.call(rbind,lapply(output, function(df) df["p.value"]))
    coefficient <- cbind(est, stdr, tvalue, p.value)
    colnames(coefficient) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    tables[[i]] <- coefficient
  }
  out <- list(call=object$call, topics=topics, tables=tables)
  class(out) <- "summary.estimateMixedEffect"
  return(out)
}

# # a function for performing mixed linear effect model
# mix.lm <- function(y, xmat, alt.formula, null.formula){
#     xdat <- xmat
#     xdat$y <- y
#     # varlist <- c(varlist, "(1|sample)")
#     alt.form <- as.formula(paste("y", 
#                                  paste(c(alt.formula,"(1|sample)"), 
#                                        collapse = "+")))
#     # random intercept 
#     int.lmer = lmer(alt.form, REML = F, data = xdat)
#     
#     # null model
#     if (!is.null(null.formula)) {
#         null.formula <- formula(paste(as.character(null.formula)[c(1,3)], collapse = " "))
#         null.form <- as.formula(paste("y", 
#                                      paste(c(null.formula,"(1|sample)"), 
#                                            collapse = "+")))
#     } else{
#         null.form <- as.formula("y ~ (1|sample)")
#     }
#     # null.x <- varlist[!(varlist %in% "time")]
#     # null.formula <- as.formula(paste("y~", 
#     #                                  paste(null.x, collapse = " + ")))
#     null.fit <- lmer(null.form, REML = F, data = xdat)
#     # perform LRT 
#     res.lrt <- anova(null.fit, int.lmer)
#     res.lm <- summary(int.lmer, ddf = "Satterthwaite")
#     out <- list(coef = res.lm$coefficients, 
#          resduals = res.lm$residuals,
#          vcov = res.lm$vcov,
#          LRT = res.lrt)
#     out
# }
    # slp.lmer = lmer(y ~ xmat$time + (1 + xmat$time|xmat$sample), REML = F)
    # testing if random slope is necessary
    # lm.res <- anova(int.lmer, slp.lmer)
    # if (lm.res$`Pr(>Chisq)`[2] < 0.05) {slope = TRUE}
    # if (slope){
    #     lmer.fit <- slp.lmer
    #     null.fit <- lmer(y ~ (1|xmat$sample), REML = F)
    #     res <- anova(null.fit, slp.lmer)
    # } else{
    #     lmer.fit <- int.lmer
    #     null.fit <- lmer(y ~ (1|xmat$sample), REML = F)
    #     res <- anova(null.fit, int.lmer)
    # }

# A function for performing simple linear regression with a cached QR decomposition
# this should be lighter weight than lm().  Works with summary.qr.lm() to give
# vcov calculations etc.
# qr.lm <- function(y, qx) {
#     if(length(y)!=nrow(qx$qr)) {
#         #probably don't match because of a prior
#         if(length(y)!=(nrow(qx$qr)-ncol(qx$qr))) stop("number of covariate observations does not match number of docs")
#         #if it passes this check its the prior. thus
#         y <- c(y,rep(0, ncol(qx$qr)))
#     }
#     beta <- solve.qr(qx, y)
#     residuals <- qr.resid(qx,y)
#     fitted.values <- qr.fitted(qx,y)
#     df.residual <- length(fitted.values) - qx$rank
#     out <- list(coefficients=beta, residuals=residuals, 
#                 fitted.values=fitted.values, 
#                 df.residual=df.residual, rank=qx$rank, qr=qx)
#     out 
# }
#this function rewrites the summary.lm() function
# to calculate from our reduced regression
# summary.qr.lm <- function (object) {
#     z <- object
#     p <- z$rank
#     rdf <- z$df.residual
#     
#     Qr <- object$qr
#     n <- nrow(Qr$qr)
#     p1 <- 1L:p
#     r <- z$residuals
#     f <- z$fitted.values
#     
#     mss <- ifelse(attr(z$terms, "intercept"), sum((f - mean(f))^2), sum(f^2)) 
#     rss <- sum(r^2)
#     
#     resvar <- rss/rdf
#     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
#     se <- sqrt(diag(R) * resvar)
#     est <- z$coefficients[Qr$pivot[p1]]
#     sigma <- sqrt(resvar)
#     list(est=est, vcov=(sigma^2 * R))
# }
