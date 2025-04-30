#Workhorse Function for the STM model
#compared to the original we have more initializations,
# more explicit options, trimmed fat, memoization

init.likelihood <- function(documents, vocab, settings, model) {
    
    globaltime <- proc.time()
    verbose <- settings$verbose
    
    ngroups <- settings$ngroups
    samples <- settings$dim$samples
    nsamples <- as.vector(table(samples))
    
    if(settings$dim$K==0) settings$dim$K <- nrow(model$beta[[1]])
    mu <- list(mu=model$mu)
    sigma <- model$sigma
    sigs <- model$sigs
    beta <- list(beta=model$beta)
    pi <- model$pi
    # beta <- model$beta
    if(!is.null(model$kappa)) beta$kappa <- model$kappa
    lambda <- model$lambda
    convergence <- NULL
    #discard the old object
    rm(model)
    
    #Pull out some book keeping elements
    ntokens <- sum(settings$dim$wcounts$x)
    betaindex <- settings$covariates$betaindex
    stopits <- FALSE
    
    # initialize sufficient statistics
    suffstats <- vector(mode="list", length=ngroups)
    
    if(settings$convergence$max.em.its==0) {
        stopits <- TRUE
        if(verbose) cat("Returning Initialization.")
    }
    # browser()
    ############
    #Step 2: Run EM
    ############
    while(!stopits) {
        t1 <- proc.time()
        #run the model
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
        pi <- suffstats$pi
        omega <- suffstats$omega
        beta.ss <- suffstats$beta
        bound.ss <- suffstats$bound
        nu <- suffstats$nu
        phi <- suffstats$phis
        #do the m-step
        mu <- opt.mu(lambda=lambda, pi = pi,
                     nsamples = nsamples, mode=settings$gamma$mode,
                     covar=settings$covariates$X, enet=settings$gamma$enet, ic.k=settings$gamma$ic.k,
                     maxits=settings$gamma$maxits)
        
        sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, omega = omega,
                           pi = pi, samples = samples,
                           mu=mu$mu, sigprior=settings$sigma$prior)
        
        beta <- opt.beta(beta.ss, beta$kappa, settings)
        sigs <- opt.sigs(pi, omega, samples)
        
        if(verbose) {
            timer <- floor((proc.time()-t1)[3])
            msg <- ifelse(timer>1,
                          sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                          "Completed M-Step. \n")
            cat(msg)
    }
        #Convergence
        # cat("Bound is ", bound.ss, "\n")
        # cat("Convergence is ", convergence, "\n")
        bound <- llh.bound(bound.ss, pi, sigs, omega, phi)
        cat("bound \n")
        convergence <- convergence.check(bound, convergence, settings)
        stopits <- convergence$stopits
        # cat("stopits is", stopits, "\n")
        
        #Print Updates if we haven't yet converged
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
    model <- list(mu=mu, sigma=sigma, beta=beta, 
                  psi = pi, sigs = sigs, settings=settings,
                  vocab=vocab, DocName = names(documents), 
                  sampleID = samples, convergence=convergence,
                  theta=exp(lambda - log(rowSums(exp(lambda)))),
                  #note altered from row.lse above because of a
                  #Windows specific bug that was happening with
                  #matrixStats package and large matrices 8/27
                  eta=lambda[,-ncol(lambda), drop=FALSE],
                  nu = nu,
                  invsigma=solve(sigma), time=time, 
                  version=utils::packageDescription("stm")$Version)
    
    class(model) <- "STM"
    return(model)
}






