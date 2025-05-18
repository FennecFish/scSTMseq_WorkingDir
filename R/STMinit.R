#Initializing Functions
#major differences here are:
# (1) better option setting
# (2) returning beta as well as mu, sigma, lambda etc.

library(fastTopics)
stm.init <- function(documents, settings) {
  
  K <- settings$dim$K
  V <- settings$dim$V
  A <- settings$dim$A
  N <- settings$dim$N
  I <- settings$dim$I
  samples <- as.factor(settings$dim$samples)
  mode <- settings$init$mode
  nits <- settings$init$nits 
  alpha <- settings$init$alpha 
  eta <- settings$init$eta 
  burnin <- settings$init$burnin 
  maxV <- settings$init$maxV
  sce <- settings$sce

  if(mode=="Spectral" | mode=="SpectralRP") {
      verbose <- settings$verbose
      if(K >= V) stop("Spectral initialization cannot be used for the overcomplete case (K greater than or equal to number of words in vocab)")
      # (1) Prep the Gram matrix
      if(verbose) cat("\t Calculating the gram matrix...\n")
      docs <- doc.to.ijv(documents) # defined in STM functions
      mat <- Matrix::sparseMatrix(docs$i,docs$j, x=docs$v) # create count matrix, row is doc, col are the vocab
      rm(docs)
      wprob <- Matrix::colSums(mat) # word count for each vocab across doc
      wprob <- wprob/sum(wprob) # prob of word for each vocab
      if(mode=="Spectral") {
          keep <- NULL
          if(!is.null(maxV)) {
              if(verbose) cat(sprintf("\t Using only %i most frequent terms during initialization...\n", maxV))
              keep <- order(wprob, decreasing=TRUE)[1:maxV]
              mat <- mat[,keep]
              wprob <- wprob[keep]
          }
          Q <- gram(mat) # defined in Spectral.R # V by V matrix
          #verify that there are no zeroes
          Qsums <- rowSums(Q)
          
          if(any(Qsums==0)) {
              #if there are zeroes, we want to remove them for just the anchor word procedure.
              temp.remove <- which(Qsums==0)
              if(is.null(keep)) {
                  keep <- which(Qsums!=0)
              } else {
                  keep <- keep[which(Qsums!=0)]
              }
              Q <- Q[-temp.remove,-temp.remove]
              Qsums <- Qsums[-temp.remove]
              wprob <- wprob[-temp.remove]
          }
          Q <- Q/Qsums
      } else {
          keep <- NULL #we have to set this because it is referenced later.
          Q <- gram.rp(mat, s=settings$init$s, p=settings$init$p,
                       d.group.size=settings$init$d.group.size, verbose=verbose)
      }
      
      
      # (2) anchor words
      if(K!=0) {
          if(verbose) cat("\t Finding anchor words...\n \t")
          anchor <- fastAnchor(Q, K=K, verbose=verbose) #spectral.R
      } else {
          if(verbose) cat("\t Finding anchor words...\n \t")
          anchor <- tsneAnchor(Q, verbose=verbose,
                               init.dims=settings$init$tSNE_init.dims,
                               perplexity=settings$init$tSNE_perplexity) #run the Lee and Mimno (2014) algorithm
          K <- length(anchor) # update K
      }
      # (3) recoverL2 #spectral.R
      if(verbose) cat("\n\t Recovering initialization...\n \t")
      # beta is a K by V matrix
      beta <- recoverL2(Q, anchor, wprob, verbose=verbose, recoverEG=settings$init$recoverEG)$A
      if(!is.null(keep)) {
          #if there were zeroes, reintroduce them
          #add in a little noise here. essentially summed over the vocab
          #it should only be .1% of the total.
          beta.new <- matrix(0, nrow=K, ncol=V)
          beta.new[,keep] <- beta
          beta.new <- beta.new + .001/V
          beta <- beta.new/rowSums(beta.new)
          rm(beta.new)
      }
      # (4) generate other parameters
      mu <- matrix(0, nrow=(K-1),ncol=1)
      sigma <- diag(5, nrow=(K-1))
      lambda <- matrix(0, nrow=N, ncol=(K-1))
      pi <- matrix(rep(rep(0, K-1), I), nrow = I)
      sigs <- diag(5, nrow = K-1)
      # omega <- diag(30, nrow=I)
      
      if(verbose) cat("Initialization complete.\n")
  }
  
  if(mode == "TopicScore") {
        cat("Initialization with topicScore. \n")
      fit <- fastTopics::fit_topic_model(t(counts(sce)),
                                          k = K,
                                         init.method = "topicscore",
                                         verbose = "none",
                                         numiter.main = 20,
                                         numiter.refine = 20)
      
        # fit <- fastTopics::init_poisson_nmf(t(counts(sce)),
        #                                    k = K, 
        #                                    init.method = "topicscore",
        #                                    verbose = "detailed",
        #                                    control = list(gc = NA))
        beta <- t(fit$F)
        beta <- beta/rowSums(beta)
        #theta <- fit$L
        # # calculate adjrandindex
        # max_indices <- apply(theta, 1, which.max)
        # fastTopics_cluster <- colnames(theta)[max_indices]
        # names(fastTopics_cluster) <- rownames(theta)
        # ft.adjr <- adjustedRandIndex(sce$Group, fastTopics_cluster[match(sce$Cell, names(fastTopics_cluster))]) 
        # cat("fastTopic adjRand is",  ft.adjr, "\n")
        # rm(theta)
        # 
        # theta <- theta/rowSums(theta)
        # lambda <- log(theta) - log(theta[,K]) #get the log-space version
        # lambda <- lambda[,-K, drop=FALSE] #drop off the last column
        # rm(theta) #clear out theta
        lambda <- matrix(0, nrow=N, ncol=(K-1))
        mu <- matrix(0, nrow=(K-1),ncol=1)
        sigma <- diag(5, nrow=(K-1))
        # mu <- colMeans(lambda) #make a globally shared mean
        # mu <- matrix(mu, ncol=1)
        # sigma <- cov(lambda)
        # pi <- rep(0,I)
        # sigs <- diag(20, nrow=I, ncol = I)
        pi <- matrix(rep(rep(0, K-1), I), nrow = I)
        sigs <- diag(5, nrow = K-1)
        # temp <- cbind(lambda, samples) %>%
        #     as.data.frame() %>%
        #     tidyr::pivot_longer(cols = !matches("^samples$"), names_to = "topic", values_to = "value") %>%
        #     group_by(samples, topic) %>%
        #     summarise(avg = mean(value), .groups = "drop") %>%
        #     tidyr::pivot_wider(names_from = topic, values_from = avg) %>%
        #     select(-samples)
        # pi <- rowMeans(temp)
        # pi <- matrix(pi, ncol = 1)
        # sigs <- apply(temp, 1, var)
        # sigs <- diag(sigs, nrow = I)
        # mu <- colMeans(lambda) #make a globally shared mean
        # mu <- matrix(mu, ncol=1)
        # sigma <- cov(lambda)  
        # rm(temp)
  }
  
  if(mode == "Random") {
      cat("Initialization with Poisson NMF \n")
      fit <- fastTopics::init_poisson_nmf(t(counts(sce)),
                                          k = K, 
                                          init.method = "random",
                                          verbose = "none")
      beta <- t(fit$F)
      beta <- beta/rowSums(beta)
      # theta <- fit$L
      # max_indices <- apply(theta, 1, which.max)
      # fastTopics_cluster <- colnames(theta)[max_indices]
      # names(fastTopics_cluster) <- rownames(theta)
      # ft.adjr <- adjustedRandIndex(sce$Group, fastTopics_cluster[match(sce$Cell, names(fastTopics_cluster))]) 
      # cat("random adjRand is",  ft.adjr, "\n")
      # rm(theta)
      # theta <- theta/rowSums(theta) # normalize theta
      # lambda <- log(theta) - log(theta[,K]) #get the log-space version
      # lambda <- lambda[,-K, drop=FALSE] #drop off the last column
      # rm(theta) #clear out theta
      
      lambda <- matrix(0, nrow=N, ncol=(K-1))
      mu <- matrix(0, nrow=(K-1),ncol=1)
      sigma <- diag(5, nrow=(K-1))
      # mu <- colMeans(lambda)
      # mu <- matrix(mu, ncol=1)
      # sigma <- cov(lambda)
      
     # patient level randomization
      # pi <- rep(0,I)
      # sigs <- diag(20, nrow=I, ncol = I)
      
      pi <- matrix(rep(rep(0, K-1), I), nrow = I)
      sigs <- diag(5, nrow = K-1)
      # temp <- cbind(lambda, samples) %>%
      #     as.data.frame() %>%
      #     tidyr::pivot_longer(cols = !matches("^samples$"), names_to = "topic", values_to = "value") %>%
      #     group_by(samples, topic) %>%
      #     summarise(avg = mean(value), .groups = "drop") %>%
      #     tidyr::pivot_wider(names_from = topic, values_from = avg) %>%
      #     select(-samples)
      # pi <- rowMeans(temp)
      # pi <- matrix(pi, ncol = 1)
      # sigs <- apply(temp, 1, var)
      # sigs <- diag(sigs, nrow = I)
      # mu <- colMeans(lambda) #make a globally shared mean
      # mu <- matrix(mu, ncol=1)
      # sigma <- cov(lambda)  
      # rm(temp)
  }

  #turn beta into a list and assign it for each aspect
  beta <- rep(list(beta),A)
  model <- list(mu=mu, sigma=sigma, sigs = sigs, 
                beta=beta, lambda=lambda, pi = pi)
  #initialize the kappa vectors
  if(!settings$kappa$LDAbeta) {
    model$kappa <- kappa.init(documents, K, V, A, interactions=settings$kappa$interactions)
  }
  
  #For custom models we already have a random init, now fill in beta
  if(mode=="Custom") {
    newbeta <- settings$init$custom
    if(!is.list(newbeta)) stop("Custom beta input is not a list")
    if(length(newbeta)!=length(beta)) stop("Custom beta list is not the same length as model specification")
    if(any(dim(newbeta[[1]])!=dim(beta[[1]]))) stop("Dimensions of custom beta do not match model specification")
    if(sum(newbeta[[1]][1,])==1) stop("It looks like you provided an unlogged version of custom beta. See documentation")
    #okay at this point we probably have checked it enough- copy it over.
    model$beta <- lapply(newbeta, exp)
  }
  
  return(model)
}

###
# Kappa initialization
###
kappa.init <- function(documents, K, V, A, interactions) {
  kappa.out <- list()
  #Calculate the baseline log-probability (m)
  freq <- matrix(unlist(documents),nrow=2) #break it into a matrix
  freq <- split(freq[2,], freq[1,]) #shift into list by word type
  m <- unlist(lapply(freq, sum)) #sum over the word types
  m <- m/sum(m)
  #m <- log(m)
  m <- log(m) - log(mean(m)) #logit of m
  kappa.out$m <- m
  
  #Defining parameters
  aspectmod <- A > 1
  if(aspectmod) {
    interact <- interactions 
  } else {
    interact <- FALSE
  }
  
  #Create the parameters object
  parLength <- K + A*aspectmod + (K*A)*interact
  kappa.out$params <- vector(mode="list",length=parLength)
  for(i in 1:length(kappa.out$params)) {
    kappa.out$params[[i]] <- rep(0, V)
  }
  
  #Create a running sum of the kappa parameters starting with m
  kappa.out$kappasum <- vector(mode="list", length=A)
  for (a in 1:A) {
    kappa.out$kappasum[[a]] <- matrix(m, nrow=K, ncol=V, byrow=TRUE)
  }
  
  #create covariates. one element per item in parameter list.
    #generation by type because its conceptually simpler
  if(!aspectmod & !interact) {
    kappa.out$covar <- list(k=1:K, a=rep(NA, parLength), type=rep(1,K))
  }
  if(aspectmod & !interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A)), a=c(rep(NA, K), 1:A), type=c(rep(1,K), rep(2,A)))      
  }
  if(interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A), rep(1:K,A)), 
                        a=c(rep(NA, K), 1:A, rep(1:A,each=K)), 
                        type=c(rep(1,K), rep(2,A), rep(3,K*A)))            
  }
  return(kappa.out)
}


