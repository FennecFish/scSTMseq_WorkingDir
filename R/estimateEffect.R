estimateEffect <- function(formula,
                           stmobj, ref.vec = NULL,# metadata=NULL, 
                           sampleNames = NULL, sampleIDs = NULL, # if id = null, then combine all samples
                           uncertainty=c("Global", "Local", "None"), documents=NULL,
                           nsims=25, prior=NULL) {
  origcall <- match.call()
  thetatype <- match.arg(uncertainty)
  # match the metadata with the docName
  metadata <- colData(stmobj$settings$sce)
  if(thetatype=="None") nsims <- 1 #override nsims for no uncertainty

  # if(!is.null(documents)) {
  #   # Convert the corpus to the internal STM format
  #   args <- asSTMCorpus(documents, data=metadata)
  #   documents <- args$documents
  #   metadata <- args$data
  # }
  
  ##
  #Step 1: Extract the formula and do some error checking
  ##
  if(!inherits(formula,"formula")) stop("formula must be a formula object.")
  if(!is.null(metadata) & !is.data.frame(metadata)) metadata <- as.data.frame(metadata)
  
  if (!is.null(ref.vec)) {
      # Extract the covariate terms from the formula
      covariate_terms <- strsplit(as.character(formula)[3], "\\s*\\*\\s*")[[1]]
      if (length(covariate_terms) != length(ref.vec)) stop("The size of the reference group should be the same as the number of covariates.")

      # Apply reference and alternative to each covariate and reference pair
      for (i in seq_along(covariate_terms)) {
        covariate <- covariate_terms[i]
        if(!covariate %in% names(ref.vec)) stop("The name of the reference group vector should be the covariate specified.")
        ref_value <- ref.vec[match(covariate, names(ref.vec))]
        if (!ref_value %in% metadata[[covariate]]) stop("The reference value '", ref_value, "' should be an element in the covariate '", covariate, "'")
        # Update the factor levels with the reference value as the first level
        metadata[[covariate]] <- factor(metadata[[covariate]], levels = c(ref_value, setdiff(unique(metadata[[covariate]]), ref_value)))
      }
  }

  termobj <- terms(formula, data=metadata)
  if(attr(termobj, "response")==1){
    #if a response is specified we have to parse it and remove it.
    # as.character of a formula turns
    # dv ~ iv
    # into:  c("~", "dv", "iv")
    response <- as.character(formula)[2] #second object is the response in this cases
    K <- eval(parse(text=response))
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
  if(!is.null(sampleIDs)) {
    if(is.null(sampleNames)) stop("Please specify the colname name for sampleIDs in the input metadata")
    if(!sampleIDs %in% stmobj$sampleID) stop("Sample ids specified must be exactly the same as the ones used in the model")
    metadata = metadata[metadata[,sampleNames] == sampleIDs,]
    subDocName <- stmobj$DocName[stmobj$sampleID %in% sampleIDs]
    stmobj <- STMsubset(stmobj, subDocName) 
  } 
  mf <- model.frame(termobj, data=metadata)
  xmat <- model.matrix(termobj,data=metadata)
  varlist <- all.vars(termobj)
  if(!is.null(metadata)) {
    data <- metadata[, varlist, drop=FALSE]
  } else {
    templist <- list()
    for(i in 1:length(varlist)) {
      templist[[i]] <- get(varlist[i])
    }
    data <- data.frame(templist)
    names(data) <- varlist
    rm(templist)
  }
  metadata <- data
  rm(data)
  ##
  #Step 2: Compute the QR decomposition
  ##
  # all the models here are essentially just OLS regressions
  # becuase we will run them many times we want to cache the 
  # expensive components in advance.
  if(!is.null(prior)) {
    if(!is.matrix(prior)) {
      prior <- diag(prior, nrow=ncol(xmat))
    } 
    if(ncol(prior)!=ncol(xmat)) stop("number of columns in prior does not match columns in design matrix")
    prior.pseudo <- chol(prior)
    xmat <- rbind(xmat,prior.pseudo)
  }
  qx <- qr(xmat)
  if(qx$rank < ncol(xmat)) {
    prior <- diag(1e-5, nrow=ncol(xmat))
    prior.pseudo <- chol(prior)
    xmat <- rbind(xmat,prior.pseudo)
    qx <- qr(xmat)
    warning("Covariate matrix is singular.  See the details of ?estimateEffect() for some common causes.
             Adding a small prior 1e-5 for numerical stability.")
  }
  ##  
  #Step 3: Calculate Coefficients
  ##
  pb <- txtProgressBar(min = 0, max = nsims, style = 3)
  storage <- vector(mode="list", length=length(K))
  for(i in 1:nsims) {
    # 3a) simulate theta
    if(thetatype=="None") thetasims <- stmobj$theta
    else {
      thetasims <- thetaPosterior(stmobj, nsims=1, type=thetatype, documents=documents)
      thetasims <- do.call(rbind, thetasims)
    }
    # logit transformation
    props.pseudo <- (thetasims + 0.01)/rowSums(thetasims + 0.01)
    thetasims <- log(props.pseudo/(1-props.pseudo))
    # 3b) calculate model
    for(k in K) {
      #lm.mod <- lm(thetasims[,k]~ xmat -1)
      #storage[[which(k==K)]][[i]] <- list(coef=coef(lm.mod),vcov=vcov(lm.mod))
      lm.mod <- qr.lm(thetasims[,k], qx)
      storage[[which(k==K)]][[i]] <- summary.qr.lm(lm.mod)      
    }
    setTxtProgressBar(pb, i)
  }

  ##
  #Step 4: Return Values
  ##
  # 4a) Give it an S3 class
  # 4b) Return the terms object, the call, the formula, the data etc.
  # 4c) anything else we need for a summary() type function
  toreturn <- list(parameters=storage, topics=K,
                   call=origcall, uncertainty=thetatype, formula=formula, data=metadata,
                   modelframe=mf, varlist=varlist, ref.vec = ref.vec)
  class(toreturn) <- "estimateEffect"
  return(toreturn)
}

# A function for performing simple linear regression with a cached QR decomposition
# this should be lighter weight than lm().  Works with summary.qr.lm() to give
# vcov calculations etc.
qr.lm <- function(y, qx) {
  if(length(y)!=nrow(qx$qr)) {
    #probably don't match because of a prior
    if(length(y)!=(nrow(qx$qr)-ncol(qx$qr))) stop("number of covariate observations does not match number of docs")
    #if it passes this check its the prior. thus
    y <- c(y,rep(0, ncol(qx$qr)))
  }
  beta <- solve.qr(qx, y)
  residuals <- qr.resid(qx,y)
  fitted.values <- qr.fitted(qx,y)
  df.residual <- length(fitted.values) - qx$rank
  out <- list(coefficients=beta, residuals=residuals, 
              fitted.values=fitted.values, 
              df.residual=df.residual, rank=qx$rank, qr=qx)
  out 
}
#this function rewrites the summary.lm() function
# to calculate from our reduced regression
summary.qr.lm <- function (object) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  
  Qr <- object$qr
  n <- nrow(Qr$qr)
  p1 <- 1L:p
  r <- z$residuals
  f <- z$fitted.values
  
  mss <- ifelse(attr(z$terms, "intercept"), sum((f - mean(f))^2), sum(f^2)) 
  rss <- sum(r^2)
  
  resvar <- rss/rdf
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  sigma <- sqrt(resvar)
  list(est=est, vcov=(sigma^2 * R))
}


#'Summary for estimateEffect
#'
#'Create a summary regression table similar to those produced for \code{lm}
#'
#'This function along with \code{\link{print.summary.estimateEffect}} creates
#'regression tables that look like typically summaries you see in R.  In general
#'we recommend that you use non-linearities such as splines via function like
#'\code{\link{s}} and in those circumstances the tables are not particularly
#'interpretable.  
#'
#'Confidence intervals are calculated by using draws from the covariance matrix
#'of each simulation to estimate the standard error.  Then a t-distribution approximation
#'is applied to calculate the various quantities of interest.
#'
#'@param object an object of class \code{"estimateEffect"}, usually a result of a call to
#'\code{\link{estimateEffect}}
#'@param topics a vector containing the topic numbers for each a summary is to be calculated.
#'Must be contained in the original \code{estimateEffect} object
#'@param nsim the number of simulations to use per parameter set to calculate the standard error.
#'Defaults to 500
#'@param ... further arguments passed to or from other methods
#'
#'@seealso \code{\link{estimateEffect}} \code{\link{plot.estimateEffect}}
#'@method summary estimateEffect  
#'@aliases summary.estimateEffect print.summary.estimateEffect
#'@export
summary.estimateEffect <- function(object, topics=NULL, nsim=500, ...) {
  if(is.null(topics)) topics <- object$topics
  if(any(!(topics %in% object$topics))) {
    stop("Some topics specified with the topics argument are not available in this estimateEffect object.")
  }
  tables <- vector(mode="list", length=length(topics))
  for(i in 1:length(topics)) {
    topic <- topics[i]
    sims <- lapply(object$parameters[[which(object$topics==topic)]], function(x) rmvnorm(nsim, x$est, x$vcov))
    sims <- do.call(rbind,sims)
    est<- colMeans(sims)
    se <- sqrt(apply(sims,2, stats::var))
    tval <- est/se
    rdf <- nrow(object$data) - length(est)
    p <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
    
    coefficients <- cbind(est, se, tval, p)
    rownames(coefficients) <- attr(object$parameters[[1]][[1]]$est, "names") 
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    tables[[i]] <- coefficients
  }
  out <- list(call=object$call, topics=topics, tables=tables)
  class(out) <- "summary.estimateEffect"
  return(out)
}

#'@method print summary.estimateEffect  
#'@export
print.summary.estimateEffect <- function(x, digits = max(3L, getOption("digits") - 3L), 
                                         signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  for(i in 1:length(x$tables)) {
    cat(sprintf("\nTopic %i:\n", x$topics[i]))
    cat("\nCoefficients:\n")
    coefs <- x$tables[[i]]
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                        na.print = "NA", ...)
    cat("\n")
  }
  invisible(x)
}