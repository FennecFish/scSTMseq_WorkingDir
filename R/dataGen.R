# This script is to generate data according to scSTMseq mode

simulate_count <- function(nSample, nCell, nGene, K, Prevalance = FALSE){
  n <- nSample*nCell
  
  # Theta is the cell by topic matrix
  # first generate sample variation
  psi <- matrix(rnorm(nSample*(K-1), mean = 0, sd = 1), nrow = nSample)
  psi <- psi[rep(1:nrow(psi), each = nCell), ]
  # then generate mu: 1) create design matrix, 2) generate gamma from a normal distribution
  # 1) Generate X, assume there exists both Response/nonResponse and Timepoints
  # Response <- rep(c(rep("nonResponse", floor(nSample / 2)), rep("Response", ceiling(nSample / 2))), each = nCell)
  Time <- rep(c("Time1", "Time2"), length.out = nCell*nSample)
  X <- data.frame(Intercept = 1,
                  Time2 = ifelse(Time == "Time2", 1, 0))
                  # Response = ifelse(Response == "Response", 1, 0))
  # X$Time2_Response = X$Time2*X$Response
  X <- as.matrix(X, ncol = ncol(X))
  # 2) generate gamma from a multivariate normal distribution
  mean_vector <- c(0.2, -0.2, rep(0, (K-3)))
  random_matrix <- matrix(rnorm((K-1)^2, mean = 0, sd = 0.05), nrow = K-1)
  cov_matrix <- crossprod(random_matrix)  # Ensure the covariance matrix is positive definite
  gamma <- mvrnorm(n = 1, mu = mean_vector, Sigma = cov_matrix)
  gamma <- rbind(rnorm(n = (K-1), mean = 0, sd = 0.1), gamma)
  # 3) generate mu
  mu = X %*% gamma
  # add mu and psi to get theta
  lambda <- mu + psi
  Theta <- exp(lambda - log(rowSums(exp(lambda))))
  # Then we want to generate beta
  
 sim_count <- matrix(as.double(rpois(n*nGene,tcrossprod(Beta, Theta))),n, nGene)
}
