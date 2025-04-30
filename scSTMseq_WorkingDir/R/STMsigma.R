#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu, omega, lambda, pi, mu, sigprior, samples) {  
    
    sigma <- list()
    I <- length(unique(samples))
    if(!is.null(pi)){
        for (i in 1:I) {
            Ni <- which(samples == unique(samples)[i])
            psi <- pi[i,]
            psi <- matrix(rep(psi, length(Ni)), nrow = length(Ni), byrow = T)
            # sigs <- diag(diag(omega)[i], nrow = nrow(mu), ncol = nrow(mu))
            #find the covariance
            if(ncol(mu)==1) {
                covariance <- crossprod(sweep(lambda[Ni,], 2, STATS=as.numeric(mu[,Ni] + psi), FUN="-"))
            } else {
                covariance <- crossprod(matrix(lambda[Ni,], nrow = length(Ni))-matrix(t(mu)[Ni,], nrow = length(Ni)) - psi) 
            }
            sigma[[i]] <- (covariance + nu)/length(Ni) #add to estimation variance
        }
        sigma <- Reduce("+", sigma)/I + omega
        sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    }else{
        #find the covariance
        if(ncol(mu)==1) {
            covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
        } else {
            covariance <- crossprod(lambda-t(mu)) 
        }
        sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
        sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
    }

    return(sigma)
}


