
stm.control <- function(documents, vocab, settings, model=NULL) {

  globaltime <- proc.time()
  verbose <- settings$verbose
  ##########
  #Step 1: Initialize Parameters
  ##########
  ngroups <- settings$ngroups
  samples <- settings$dim$samples
  nsamples <- as.vector(table(samples))
  I <- settings$dim$I
  if(is.null(model)) {
    if(verbose) cat(switch(EXPR=settings$init$mode,
                           Spectral = "Beginning Spectral Initialization \n",
                           LDA = "Beginning LDA Initialization \n",
                           Random = "Beginning Random Initialization \n",
                           Custom = "Beginning Custom Initialization \n"))
    #initialize
    
    model <- stm.init(documents, settings) # see STMinit.R
    #if we were using the Lee and Mimno method of setting K, update the settings
    if(settings$dim$K==0) settings$dim$K <- nrow(model$beta[[1]])
    #unpack
    
    mu <- list(mu=model$mu)
    sigma <- model$sigma
    beta <- list(beta=model$beta)
    # if only one sample, makes sample variation null
    if(I == 1) {
        pi <- NULL
        sigs <- NULL
    } else{
        pi <- model$pi
        sigs <- model$sigs
    }
    # beta <- model$beta
    if(!is.null(model$kappa)) beta$kappa <- model$kappa
    lambda <- model$lambda
    convergence <- NULL
    #discard the old object
    rm(model)
  } else {
    if(verbose) cat("Restarting Model...\n")
    #extract from a standard STM object so we can simply continue.
    mu <- model$mu
    beta <- list(beta=lapply(model$beta$logbeta, exp))
    if(!is.null(model$beta$kappa)) beta$kappa <- model$beta$kappa
    sigma <- model$sigma
    lambda <- model$eta
    # if only one sample, makes sample variation null
    if(I == 1) {
        pi <- NULL
        sigs <- NULL
    } else{
        pi <- model$pi
        sigs <- model$sigs
    }
    convergence <- model$convergence
    #manually declare the model not converged or it will stop after the first iteration
    convergence$stopits <- FALSE
    convergence$converged <- FALSE
    #iterate by 1 as that would have happened otherwise
    convergence$its <- convergence$its + 1
  }

  #Pull out some book keeping elements
  ntokens <- sum(settings$dim$wcounts$x)
  betaindex <- settings$covariates$betaindex
  stopits <- FALSE
  if(ngroups!=1) {
    # randomly assign groups so that subsample are representative
    groups <- base::split(seq_len(length(documents)),
                    sample(rep(seq_len(ngroups), length=length(documents))))
  }
  suffstats <- vector(mode="list", length=ngroups)
  # for ARI
  ARI <- vector(mode = "list")
  if(settings$convergence$max.em.its==0) {
    stopits <- TRUE
    if(verbose) cat("Returning Initialization.")
  }
 
  ############
  #Step 2: Run EM
  ############
  while(!stopits) {

      t1 <- proc.time()
      #run the model
      # source("R/STMestep.R")

      suffstats <- estep(documents=documents, beta.index=betaindex,
                         update.mu=(!is.null(mu$gamma)),
                         beta = beta, lambda.old = lambda,
                         mu = mu$mu, sigma = sigma,
                         sigs = sigs, pi.old = pi,
                         samples = samples, verbose)
      msg <- sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
      if(verbose) cat(msg)
      t1 <- proc.time()
      sigma.ss <- suffstats$sigma
      lambda <- suffstats$lambda
      beta.ss <- suffstats$beta
      bound.ss <- suffstats$bound
      nu <- suffstats$nu
      phi <- suffstats$phis
      if(!is.null(pi)){
          pi <- suffstats$pi
          omega <- suffstats$omega
      }
     
      # trace.ss <- suffstats$trace
      # new_bound.ss <- suffstats$new_bound
      
      #do the m-step

      mu <- opt.mu(lambda=lambda, pi = pi,
                   nsamples = nsamples, mode=settings$gamma$mode,
                   covar=settings$covariates$X, enet=settings$gamma$enet, ic.k=settings$gamma$ic.k,
                   maxits=settings$gamma$maxits, settings = settings)
      beta <- opt.beta(beta.ss, beta$kappa, settings)

      # sce_aggregated <- aggregateAcrossCells(sce, ids = colData(sce)$Batch)
      # set <- newSeqExpressionSet(counts(sce_aggregated))
      # y <- DGEList(counts=counts(set))
      # y <- calcNormFactors(y, method="upperquartile")
      # y <- estimateGLMCommonDisp(y)
      # y <- estimateGLMTagwiseDisp(y)
      # 
      # fit <- glmFit(y)
      # res <- residuals(fit, type="deviance")
      # set4 <- RUVr(set, rownames(sce_aggregated), k=1, res)
      # pi <- matrix(rep(pData(set4)[,1], K-1), nrow = I)

      
      if(!is.null(pi)){
          sigs <- opt.sigs(pi, omega, samples)
          alpha <- pi
          sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, omega = omega,
                             pi = alpha, samples = samples,
                             mu=mu$mu, sigprior=settings$sigma$prior)
      } else{

          sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, omega = NULL,
                             pi = NULL, samples = samples,
                             mu=mu$mu, sigprior=settings$sigma$prior)
      }
      timer <- floor((proc.time()-t1)[3])
      msg <- sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
      if(verbose) cat(msg)
    #Convergence
    # cat("Bound is ", bound.ss, "\n")
    # cat("Convergence is ", convergence, "\n")

      # this part is for testing purpose only (calculating ARI)
      lambda_temp <- cbind(lambda,0)
      theta_temp <- exp(lambda_temp - log(rowSums(exp(lambda_temp))))
      # settings$sce$Group
      max_indices <- apply(theta_temp, 1, which.max)
      colnames(theta_temp) <- paste0("topic_", 1:ncol(theta_temp))
      res_cluster <- colnames(theta_temp)[max_indices]
      names(res_cluster) <- names(documents)
      true_group <- settings$sce$Group[match(names(res_cluster), settings$sce$Cell)]
      result <- tryCatch({
        ARI <- append(ARI, adjustedRandIndex(res_cluster, true_group))
        TRUE  # No error occurred
      }, error = function(e) {
        cat("The ARI has failed with ARI =", ARI, "\n")
        cat("True Group is ", true_group, "\n")
        cat("Inferred Cluster is ", res_cluster, "\n")
        # cat("Appended Item is ", adjustedRandIndex(res_cluster, true_group), "\n")
        # browser()  # Enters debug mode
        FALSE  # Indicates an error occurred
      })
      # ARI <- append(ARI, adjustedRandIndex(res_cluster, true_group))
      
      # browser()
    # bound <- llh.bound(bound.ss, alpha, sigs, omega, phi)
    bound <- sum(bound.ss)
    # trace <- sum(trace.ss)
    # new_bound <- sum(new_bound.ss)
    # 
    cat("calculate log likelihood \n")
    # cat("bound \n")
    convergence <- convergence.check(bound, convergence, settings)
    stopits <- convergence$stopits
    # cat("stopits is", stopits, "\n")
    #Print Updates if we haven't yet converged
    # if(!stopits & verbose) report(convergence, ntokens=ntokens, beta, vocab,
    #                                    settings$topicreportevery, verbose)
    if(!stopits & verbose) report(convergence, ntokens=ntokens, beta, vocab,
                                  settings$topicreportevery, verbose)
  }
  #######
  #Step 3: Construct Output
  #######
  time <- (proc.time() - globaltime)[3]
  #convert the beta back to log-space
  beta$logbeta <- beta$beta
  for(i in 1:length(beta$logbeta)) {
    beta$logbeta[[i]] <- safelog(beta$logbeta[[i]])
  }
  beta$beta <- NULL
  lambda <- cbind(lambda,0)
  if(is.null(pi)){
      alpha <- NULL
      sigs <- NULL
  }

  model <- list(mu=mu, sigma=sigma, beta=beta, 
                psi = list(alpha = alpha, sigs = sigs), settings=settings,
                vocab=vocab, DocName = names(documents), 
                sampleID = samples, convergence=convergence,
                theta=exp(lambda - log(rowSums(exp(lambda)))),
                eta=lambda[,-ncol(lambda), drop=FALSE],
                nu = nu,
                ARI = ARI,
                time=time, 
                version=utils::packageDescription("stm")$Version)

  class(model) <- "STM"
  return(model)
}






