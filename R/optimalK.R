
optimal_K <- function(sce, K, init.type = "Spectral", sample = NULL,
                    N=floor(.25*ncol(sce)), proportion = .5, 
                    heldout.seed = NULL, M = 10, cores = 1, ...) {

  #Make a heldout dataset
  heldout <- make.heldout.sce(sce, N=N, proportion=proportion, 
                          seed=heldout.seed)
  cat("Heldout Set Created \n")
  # warnings
  # if( "content" %in% names(list(...)) ) {
  #   warning("Exclusivity calculation only designed for models without content covariates", call.=FALSE)
  # }
  
  # single core
  if (cores == 1) {
      g <- list()
      for (i in seq_along(K)) { # loop produces nicer printout than lapply
          g[[i]] <- get_statistics(K[i], heldout=heldout, init.type=init.type,M=M,...)
      }
  # multi core
  } else {
      cat("Using multiple-cores.  Progress will not be shown. \n")
      g <- parallel::mclapply(K, get_statistics, mc.cores = cores, heldout=heldout, init.type=init.type,
                              M=M,...)
  } 

  # output
  res <- as.data.frame(do.call('rbind', g)) 
  res <- do.call(cbind, lapply(res, function(x) do.call(rbind, x))) %>%
      as.data.frame()
  colnames(res) <- c("K","perplexity", "heldout",
                             "residual", "bound", "lbound", 
                             "em.its")
  
  # calcualte rpc
  diff_k <- diff(res$K)
  diff_perplexity <- diff(res$perplexity)
  rpc <- abs(diff_perplexity / diff_k)
  
  res$rpc <- c(NA, rpc)
  
  toreturn <- list(results=res, call=match.call(expand.dots=TRUE))
  class(toreturn)<- "searchK"
  return(toreturn)
}

# compute statistics for a particular number of topics k
get_statistics <- function(k, heldout, init.type, M, ...) { # k = one particular topic number; K = vector of topic numbers
  out <- NULL # output vector
  out[['K']] <- k
  #run scstm
  sub_sce <- heldout$training.sce
  # args <- prepsce(sub_sce)
  # documents <- args$documents
  # vocab <- args$vocab
  # data <- args$meta
  # sub_sce <- args$sce
  model <- multi_stm(sce = sub_sce, # documents = documents, vocab = vocab, data = data,
               K=k, init.type=init.type, heldout = TRUE, ...)
  # #calculate values to return
  # if( !"content" %in% names(list(...)) ) {  # only calculate exclusivity for models without content covariates
  #   out[['exclus']] <- mean(unlist(exclusivity(model, M=M, frexw=.7)))
  #   out[['semcoh']] <- mean(unlist(semanticCoherence(model, heldout$documents, M)))
  # } 
  eval.res <- eval.heldout(model, heldout$heldout.dat)
  out[['perplexity']] <- eval.res$expected.perplexity
  out[['heldout']] <- eval.res$expected.heldout
  out[['residual']] <- checkResiduals(model,heldout)$dispersion
  out[['bound']] <- max(model$convergence$bound)
  out[['lbound']] <- max(model$convergence$bound) + lfactorial(model$settings$dim$K)
  out[['em.its']] <- length(model$convergence$bound)
  cat("Complete Topic", k, "\n")
  return(out)
}

eval.heldout <- function(model, missing) {
    perplexity <- vector(length=length(missing$index))
    heldout.likelihood <- vector(length=length(missing$index))
    # ntokens <- vector(length=length(missing$index))
    beta <- lapply(model$beta$logbeta, exp)
    bindex <- model$settings$covariates$betaindex[missing$index]
    for(i in 1:length(missing$index)) {
        cell_id <- missing$index[i]
        genes <- names(missing$cells[[i]])
        probs <- model$theta[cell_id,]%*%beta[[bindex[i]]][,match(genes, model$vocab)]
        probs <- rep(probs, missing$cells[[i]])
        heldout.likelihood[i] <- mean(log(probs))
        perplexity[i] <- exp(-mean(log(probs)))
    }
    out <- list(expected.perplexity=mean(perplexity, na.rm=TRUE), 
                expected.heldout=mean(heldout.likelihood, na.rm=TRUE)) #the mean here to deal with 0 length docs
    return(out)
}

