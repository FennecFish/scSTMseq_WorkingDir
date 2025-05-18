# calculate bound for convergence
llh.bound <- function(leta, alpha, sigs, omega, phi) {
  browser()
    leta <- sum(leta)
    # compute chol of sigma^s
    det_sigs <- det(sigs)
    # compute inv of sigma^s
    inv_sigs <- diag(1 /diag(sigs), nrow = nrow(sigs))
    tr <- sum(diag(omega %*% inv_sigs))
    Ep_psi <- -0.5*det_sigs -0.5*t(alpha) %*% inv_sigs %*% alpha - 0.5*tr
    Eq_psi <- -0.5*sum(diag(omega))
    # adding 0.1 to avoid 0
    logphi <- log(phi+ 0.1)
    Eq_z <- sum((phi+ 0.1) * logphi)
    bound <- leta + Ep_psi - Eq_psi - Eq_z
    # bound <- leta + Ep_psi - Eq_psi 
    return(bound)
}