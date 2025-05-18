checkResiduals <- function(stmobj, heldout, tol=.01) {
    missing.index <- heldout$heldout.dat$missing.gene.index
    cell <- heldout$heldout.dat$cells
    beta <- lapply(stmobj$beta$logbeta, exp)
    theta <- stmobj$theta
    index <- stmobj$settings$covariates$betaindex
    K <- stmobj$settings$dim$K
    n <- length(heldout$heldout.dat$cells)
    phat <- stmobj$settings$dim$V
    d <- n*(K-1) + K*( phat-1 )
    
    doc.resid <- function(cell, index, theta, beta) {
        q <- theta%*%beta 
        m <- sum(cell)
        Nhat <- sum(q*m > tol)
        x <- rep(0, ncol(beta))
        x[index] <- as.vector(cell)
        out <- sum((x^2 - 2*x*q*m)/(m*q*(1-q))) + sum(m*q/(1-q))
        return(list(out=out, Nhat=Nhat))
    }
    result <- vector(mode="list", length=length(cell))
    Nhat <- 0
    for(i in 1:length(result)) {
        resid <- doc.resid(cell[[i]], missing.index[[i]], theta[i,], beta[[index[i]]])
        result[[i]] <- resid$out
        Nhat <- Nhat + resid$Nhat
    }
    D <- sum(unlist(result))
    df <- Nhat - phat  - d
    sig2 <- D/df
    
    rho <- suppressWarnings(pchisq(D, df=df, lower.tail=FALSE))
    D <- list(dispersion=sig2, pvalue=rho, df=df)
    return(D) 
}
