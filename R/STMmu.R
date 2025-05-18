## Optimization for Global Parameters over Doc-Topic Proportions
#Variational Linear Regression with a Half-Cauchy hyperprior 
# (Implementation based off the various LMM examples from Matt Wand)
# This code is intended to be passed a Matrix object
# vb.variational.reg <- function(Y,X, b0=1, d0=1, Xcorr=NULL, maxits=1000) {
#     if(is.null(Xcorr)) Xcorr <- crossprod(X)
#     XYcorr <- crossprod(X,Y) 
#     
#     an <- (1 + nrow(X))/2
#     D <- ncol(X)
#     N <- nrow(X)
#     w <- rep(0, ncol(X))
#     error.prec <- 1 #expectation of the error precision
#     converge <- 1000
#     cn <- ncol(X) # - 1 for the intercept and +1 in the update cancel
#     dn <- 1
#     Ea <- cn/dn #expectation of the precision on the weights
#     ba <- 1
#     
#     ct <- 1
#     while(converge>.0001) {
#         w.old <- w
#         
#         #add the coefficient prior.  Form depends on whether X is a Matrix object or a regular matrix.
#         if(is.matrix(X)) {
#             ppmat <- diag(x=c(0, rep(as.numeric(Ea), (D-1))),nrow=D) 
#         } else {
#             ppmat <- Diagonal(n=D, x=c(0, rep(as.numeric(Ea), (D-1))))
#         }
#         invV <- error.prec*Xcorr + ppmat
#         #if its a plain matrix its faster to use the cholesky, otherwise just use solve
#         if(is.matrix(invV)) {
#             V <- chol2inv(chol(invV))
#         } else {
#             #Matrix package makes this faster even when its non-sparse
#             V <- solve(invV)     
#         }
#         w <- error.prec*V%*%XYcorr
#         
#         # parameters of noise model (an remains constant)
#         sse <- sum((X %*% w - Y)^ 2)
#         bn <- .5*(sse + sum(diag(Xcorr%*%V))) + ba
#         error.prec <- an/bn
#         ba <- 1/(error.prec + b0)
#         
#         #subtract off the intercept while working out the hyperparameters
#         # for the coefficients
#         w0 <- w[1]
#         w <- w[-1]
#         da <- 2/(Ea + d0)
#         dn <- 2*da + (crossprod(w) + sum(diag(V)[-1]))
#         Ea <- cn / dn
#         #now combine the intercept back in 
#         w <- c(w0,w)
#         ct <- ct + 1
#         if(ct > maxits) {
#             stop("Prevalence regression failing to converge within iteration limit.  May want to try gamma.prior='L1'. You can change max iterations using control.  See stm documentation")
#         }
#         converge <- sum(abs(w-w.old))
#     }
#     return(list(w = w))
# }

library(lme4)
vb.variational.reg <- function(Y,X, b0=1, d0=1, Xcorr=NULL, maxits=1000) {
    if(is.null(Xcorr)) Xcorr <- crossprod(X)
    XYcorr <- crossprod(X,Y)

    N <- nrow(X)
    D <- ncol(X)
    an <- 1 + N/2 #an is changed slightly from the original STM code
    w <- rep(0, ncol(X))
    # error.prec <- 1 #expectation of the error precision
    converge <- 1000
    cn <- 1 + D/2 # this is also changed
    dn <- 1
    Ea <- cn/dn #expectation of the precision on the weights
    b0 <- 1

    ct <- 1

    while(converge>.0001) {

        w.old <- w

        #add the coefficient prior.  Form depends on whether X is a Matrix object or a regular matrix.
        if(is.matrix(X)) {
            ppmat <- diag(x=c(0, rep(as.numeric(Ea), (D-1))),nrow=D)
        } else {
            ppmat <- Diagonal(n=D, x=c(0, rep(as.numeric(Ea), (D-1))))
        }
        invV <- Xcorr + ppmat # removed error.prec
        #if its a plain matrix its faster to use the cholesky, otherwise just use solve
        if(is.matrix(invV)) {
            V <- chol2inv(chol(invV))
        } else {
            #Matrix package makes this faster even when its non-sparse
            V <- solve(invV)
        }
        # w <- error.prec*V%*%XYcorr
        w <- V%*%XYcorr
        # parameters of noise model (an remains constant)
        sse <- sum((X %*% w - Y)^ 2)       
        # bn <- .5*(sse + sum(diag(Xcorr%*%V))) + b0
        bn <- .5*(sse + crossprod(w)*Ea) + b0
        # error.prec <- an/bn
        #subtract off the intercept while working out the hyperparameters
        # for the coefficients
        w0 <- w[1]
        w <- w[-1]
        
        eTauWTW2 <- sum(diag(V)[-1]) + crossprod(w)* an/bn
        dn <- d0 + eTauWTW2/2
        # ba <- 1/(error.prec + b0)
        
        #da <- 2/(Ea + d0)
        #dn <- 2*da + (crossprod(w) + sum(diag(V)[-1]))
        Ea <- cn / dn
        # sn <- (0.5 + sum(w^2))/(length(w) + 0.5)

        #now combine the intercept back in
        w <- c(w0,w)
        ct <- ct + 1
        if(ct > maxits) {
            stop("Prevalence regression failing to converge within iteration limit.  May want to try gamma.prior='L1'. You can change max iterations using control.  See stm documentation")
        }
        converge <- sum(abs(w-w.old))
    }
    error.prec <- as.numeric(an/bn)
    sn <- V/error.prec
    return(list(w = w, invV = invV, V = V, error.prec = error.prec, an = an))
}

#main method up top, regression-implementations below.
opt.mu <- function(lambda, pi, nsamples,
                   mode=c("Pooled", "L1", "LinearRegression", "LinearMixed"), covar=NULL, enet=NULL, ic.k=2,
                   maxits=1000, settings = NULL) {
  # #When there are no covariates we use the CTM method
  # if(mode=="CTM") {
  #   mu <- matrix(rowMeans(lambda), ncol=1)
  #   return(list(mu=mu, gamma=NULL))
  # }
  #Variational Linear Regression with a Gamma hyperprior
  if(mode=="Pooled") {
    gamma <- vector(mode="list",length=ncol(lambda))
    param <- vector(mode="list",length=ncol(lambda))
    Xcorr <- crossprod(covar)
    if(!is.null(pi)){
        pis <- pi[rep(1:nrow(pi), times = nsamples), ]
        Y_all <- lambda - pis} else{Y_all <- lambda}
    for (i in 1:ncol(lambda)) {
      # vb.res <- vb.variational.reg(Y=lambda[,i]-pi[,i], X=covar, Xcorr=Xcorr, maxits=maxits) 
        vb.res <- vb.variational.reg(Y = Y_all[,i], X=covar, Xcorr=Xcorr, maxits=maxits)
      gamma[[i]] <- vb.res$w
      param[[i]] <- vb.res
    }
    gamma <- do.call(cbind,gamma)
    rownames(gamma) <- colnames(covar)
    # sn <- do.call(cbind,sn)
    mu<- t(covar%*%gamma)
    #if its not a regular matrix,coerce it as it won't be sparse.
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma, param = param))
  }
  
  if(mode == "LinearRegression"){
    if(!is.null(pi)){
      pis <- pi[rep(1:nrow(pi), times = nsamples), ]
      Y_all <- lambda - pis
    }else{Y_all <- lambda}
    X <- covar[,-1]
    models <- vector(mode = "list")
    std.models <- vector(mode = "list")
    lm.model <- vector(mode = "list")
    for (k in seq_len(ncol(Y_all))) {
      lm.model.temp <- lm(Y_all[, k] ~ X )
      models[[k]] <-  summary(lm.model.temp)$coefficients[,"Estimate"]
      std.models[[k]] <- summary(lm.model.temp)$coefficients[,"Std. Error"]
      lm.model[[k]] <- lm.model.temp
    }
    gamma <- do.call(cbind, models)
    std.gamma <- do.call(cbind, std.models)
    mu <- t(covar%*%gamma)

    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma, std.gamma = std.gamma, lm.model = lm.model))
  }
  
  if(mode == "LinearMixed"){
    
    if(!is.null(pi)){
      pis <- pi[rep(1:nrow(pi), times = nsamples), ]
      Y_all <- lambda - pis
    }else{Y_all <- lambda}

    # check for random effect
    terms <- unlist(strsplit(deparse(settings$covariates$formula), "\\+")) 
    random_effects <- gsub("\\s+", "", terms[grepl("\\|", terms)])
    random_effects <- gsub(".*\\|(.*)\\)", "\\1", random_effects)
    
    # match Sample from sce to covar

    mixed.covar <- covar %>%
      as.data.frame() %>%
      mutate(Sample = settings$sce$Sample[match(rownames(covar), settings$sce$Cell)])
    X <- mixed.covar #[,-1]
    models <- vector(mode = "list")
    mu <- vector(mode = "list")
    lm.model <- vector(mode = "list")
    lmer.formula <- paste(c(colnames(covar)[-1], "(1|Sample)"), collapse = "+")
    lmer.formula <- as.formula(paste0("Y_all[, k]~", lmer.formula))
    for (k in seq_len(ncol(Y_all))) {
      # lm.model.temp <- lme4::lmer(Y_all[, k] ~ TimeTime2 + (1|Sample), data = X)
      lm.model.temp <- lme4::lmer(lmer.formula, data = X)
      # if(isSingular(lm.model.temp)){
      #   lm.model.temp <- lm(Y_all[, k] ~ TimeTime2, data = X)
      # }
      models[[k]] <-  summary(lm.model.temp)$coefficients[,"Estimate"]
      mu[[k]] <-predict(lm.model.temp)
      # std.models[[k]] <- summary(lm.model.temp)$coefficients[,"Std. Error"]
      lm.model[[k]] <- lm.model.temp
    }
    gamma <- do.call(cbind, models)
    mu <- t(do.call(cbind,mu))
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma, lm.model = lm.model))
  }
  
  if(mode=="L1") {

    if(!is.null(pi)){
      pis <- pi[rep(1:nrow(pi), times = nsamples), ]
      Y_all <- lambda - pis
    }else{Y_all <- lambda}
    X <- covar[,-1]
    # split into train and test 70-30
    train.index <- sample(1:nrow(Y_all), round(0.7 * nrow(Y_all)), replace = F)
    # test.index <- setdiff(1:nrow(Y_all), train.index)
    
    Y.train <- Y_all[train.index,]
    X.train <- X[train.index,]
    Y.test <- Y_all[-train.index,]
    X.test <- X[-train.index,]

    # try to find the optimal enet
    models <- list()
    for (i in 0:5) {
      name <- paste0("alpha", i/5)
      
      models[[name]] <-
        glmnet::cv.glmnet(X.train, Y.train, type.measure="mse", alpha=i/5, 
                  family="mgaussian")
    }

    results <- data.frame()
    for (i in 0:5) {
      name <- paste0("alpha", i/5)
      ## Use each model to predict 'y' given the Testing dataset
      predicted <- predict(models[[name]], s=models[[name]]$lambda.min, newx=X.test)
      predicted <- predicted[,,1]
      ## Calculate the Mean Squared Error...
      mse <- mean((Y.test - predicted)^2)
      
      ## Store the results
      temp <- data.frame(alpha=i/5, mse=mse, name=name)
      results <- rbind(results, temp)
    }
    
    # find the alpha with smallest alpha
    alpha = results$alpha[order(results$mse, decreasing = F)][1]
    out <- glmnet::glmnet(x=X, y = Y_all, family="mgaussian", 
                          alpha=alpha, lambda = models[[paste0("alpha", alpha)]]$lambda.min)

    unpack <- unpack.glmnet(out, ic.k=ic.k)
    gamma <- rbind(unpack$intercept, unpack$coef)
    mu <- t(covar%*%gamma)
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma))
  }
}



