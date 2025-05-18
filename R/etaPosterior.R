# In this script, we create a function to sample from
# the posterior distribution of Eta
# which is equvivalent to sampling from the variaitonal distribution of eta
etaPosterior <- function(model, nsims = 100){
  lambda <- model$eta
  nu <- model$nu
  
  out <- vector("list", nrow(lambda))
  row_names <- rownames(lambda)
  for (i in seq_len(nrow(lambda))) {
    # sigma <- diag(diag(nu[[i]]))
    sigma <- nu[[i]]
    lambda_i <- lambda[i,]
    out[[i]] <- mvrnorm(nsims, mu = lambda_i, Sigma = sigma)
  }
  
  SimEta <- vector("list", nsims)
  for (i in seq_len(nsims)) {
    # Create a matrix by extracting the i-th row from each list element in 'out'
    SimEta[[i]] <- matrix(unlist(lapply(out, function(mat) mat[i, ])), 
                          nrow = length(out), byrow = TRUE)
    rownames(SimEta[[i]]) <- row_names
  }
  return(SimEta)
}
