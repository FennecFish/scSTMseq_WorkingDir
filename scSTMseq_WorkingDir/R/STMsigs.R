#Optimization for the patient hereogeneity Covariance Matrix
opt.sigs <- function(pi, omega, samples) {  
    sigs <- matrix(0, nrow = nrow(omega), ncol = ncol(omega))
    I <- length(unique(samples))
    # psi <- numeric(length(unique(samples)))
    # psi <- do.call(rbind, pi)
    for (i in 1:I) {
        # Ni <- which(samples == unique(samples)[i])
        # psi[i] <- mean(colMeans(pi[Ni,]))
        sigs <- sigs + tcrossprod(pi[i,]) + omega
    }
    sigs <- sigs/I
    return(sigs)
}


