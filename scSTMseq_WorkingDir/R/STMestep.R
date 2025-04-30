#E-Step for a Document Block
#[a relatively straightforward rewrite of previous
# code with a focus on avoiding unnecessary computation.]

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                       beta, lambda.old, mu, sigma,
                       samples, pi.old, sigs,
                       verbose) {
  
  #quickly define useful constants
  V <- ncol(beta$beta[[1]])
  K <- nrow(beta$beta[[1]])
  N <- length(documents)
  A <- length(beta$beta)
  I <- length(unique(samples))
  ctevery <- ifelse(N>100, floor(N/100), 1)

  if(!update.mu) mu.l <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  
  # pi.ss <- numeric(N)
  pi.ss <- vector("list", length=N)
  pi <- vector("list", length=I)
  bound <- vector(length=N)
  #trace <- vector(length=N)
  #new_bound <- vector(length=N)
  lambda <- vector("list", length=N)
  nu <- vector(mode="list", length=N)
  for(i in 1:N) {
      nu[[i]] <- matrix(0, nrow=K-1,ncol=K-1)
  }
  # 2) Precalculate common components
    # calculate inverse and entropy of sigmat
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(inherits(sigobj,"try-error")) {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
  }
  
  # calculate inverse of sigma^s when sigma^s is not null
  if(I != 1){
      sigs_obj <- try(chol.default(sigs), silent=TRUE)
      if(inherits(sigs_obj,"try-error")) {
          sigsentropy <- (.5*determinant(sigs, logarithm=TRUE)$modulus[1])
          sigs_inv <- solve(sigs)
      } else {
          sigsentropy <- sum(log(diag(sigs_obj)))
          sigs_inv <- chol2inv(sigs_obj)
      }
      
      omega <- matrix(0, nrow = I, ncol = I)
      # calcualte omega
      sum_sig <- siginv + sigs_inv
      sum_sig_obj <- try(chol.default(sum_sig), silent=TRUE)
      if(inherits(sum_sig_obj,"try-error")) {
          omega <- solve(sum_sig)
      } else {
          omega <- chol2inv(sum_sig_obj)
      }
      
      omegaobj <- try(chol.default(omega), silent=TRUE)
      if(inherits(omegaobj,"try-error")) {
          omegaentropy <- (.5*determinant(omega, logarithm=TRUE)$modulus[1])
      } else {
          omegaentropy <- sum(log(diag(omegaobj)))
      }
      
      for (i in 1:I) {
          # psi.i <- rep(pi.old[i], ncol(lambda.old)) # repeat pi into a K-1 dimensional vector
          psi.i <- pi.old[i,]
          # sigs.i <- diag(sigs)[i]
          
          # omega <- solve(siginv + sigs_inv)
          Ni <- which(samples == unique(samples)[i])
          
          for(l in Ni) {
              # print(l)
              #update components
              doc <- documents[[l]]
              words <- doc[1,]
              aspect <- beta.index[l]
              init <- lambda.old[l,]
              if(update.mu) mu.l <- mu[,l]
              beta.l <- beta$beta[[aspect]][,words,drop=FALSE]
              
              doc.results <- logisticnormalcpp(eta=init, mu=mu.l, 
                                               psi = psi.i, 
                                               siginv=siginv, sigmaentropy=sigmaentropy,
                                               sigs = sigs, sigs_inv = sigs_inv,
                                               sigsentropy = sigsentropy,
                                               beta=beta.l, doc=doc,
                                               omega = omega, omegaentropy = omegaentropy)

              if(any(is.na(doc.results))){
                sigma.ss <- sigma.ss 
                bound[l] <- 0
                nu[[l]] <- doc.results$eta$nu
              } else{
                # update sufficient statistics 
                sigma.ss <- sigma.ss + doc.results$eta$nu
                beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
                phis <- doc.results$phis
                lambda[[l]] <- c(doc.results$eta$lambda)
                bound[l] <- doc.results$bound
                nu[[l]] <- doc.results$eta$nu
                pi.ss[[l]] <- c(doc.results$pi)
              }
       
              #trace[l] <- doc.results$trace
              #new_bound[l] <- doc.results$new_bound
              
              # adding 0.01 to avoid 0
              # logphi <- log(phis+ 0.01)
              # Eq_z <- sum(phis * logphi)
              
              #  det_sigs <- det(sigs)
              #  inv_sigs <- diag(1 /diag(sigs), nrow = nrow(sigs))
              #  tr <- sum(diag(omega %*% inv_sigs))
              #  Ep_psi <- -0.5*det_sigs -0.5*t(alpha) %*% inv_sigs %*% alpha - 0.5*tr
              #  Eq_psi <- -0.5*sum(diag(omega))
              # # bound[l] <- doc.results$bound + Ep_psi - Eq_psi
              #  bound[l] <- Ep_psi - Eq_psi
              if(verbose && l%%ctevery==0) cat(".")
          }

          pi.i <- pi.ss[Ni]
          pi[[i]] <- colMeans(do.call(rbind, pi.i)) 
          # omega[i,i] <- omega.i
      }
      # 3) Document Scheduling
      # For right now we are just doing everything in serial.
      # the challenge with multicore is efficient scheduling while
      # maintaining a small dimension for the sufficient statistics.
      
      
      
      if(verbose) cat("\n") #add a line break for the next message.
      lambda <- do.call(rbind, lambda)
      pi <- do.call(rbind, pi)
      return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, 
                  lambda=lambda, nu = nu, pi = pi, omega = omega, phis = phis))
  } else {
      ##### this section is for a single sample
      for(l in seq_along(samples)) {
          # print(l)
          #update components
          doc <- documents[[l]]
          words <- doc[1,]
          aspect <- beta.index[l]
          init <- lambda.old[l,]
          if(update.mu) mu.l <- mu[,l]
          beta.l <- beta$beta[[aspect]][,words,drop=FALSE]
          
          doc.results <- logisticnormalcpp(eta=init, mu=mu.l, 
                                           psi = NULL, omega = NULL,
                                           siginv=siginv, sigs = NULL,
                                           beta=beta.l,
                                           doc=doc,
                                           sigmaentropy=sigmaentropy)
          # update sufficient statistics 
          
          sigma.ss <- sigma.ss + doc.results$eta$nu
          beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
          phis <- doc.results$phis
          lambda[[l]] <- c(doc.results$eta$lambda)
          names(lambda)[l] <- names(documents)[l]
          bound[l] <- doc.results$bound
          nu[[l]] <- doc.results$eta$nu
          if(verbose && l%%ctevery==0) cat(".")
      }
  }
  # 3) Document Scheduling for single patient
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.

  if(verbose) cat("\n") #add a line break for the next message.
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, 
              lambda=lambda, nu = nu, phis = phis))
}

