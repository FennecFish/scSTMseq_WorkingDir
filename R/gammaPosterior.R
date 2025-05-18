# This script sample from gamma posterior distribution
# and map the gamma back to simplex
# calculate relative difference proportion change for each topic
library(mvtnorm)

gammaPosterior <- function(model, nsims = 100){

  if(model$settings$gamma$mode == "Pooled"){
    gamma_mean <- model$mu$gamma[-1,]
    invV <- lapply(model$mu$param, `[[`, "invV")[[1]][-1,-1]
    error.prec <- do.call(c, lapply(model$mu$param, `[[`, "error.prec"))
    df <- lapply(model$mu$param, `[[`, "an")[[1]]
    # gamma_shape <- 
    # gamma_variance <- model$mu$sn
    # gamma_variance <- lapply(gamma_variance, function(mat) mat[2:nrow(mat), 2:ncol(mat)])
    if(dim(model$settings$covariates$X)[2] > 2){
      sim_value <- lapply(1:ncol(gamma_mean), function(i) {
        invShape <- invV*error.prec[i]
        shape <- solve(invShape)
        x <- rmvt(nsims, sigma=shape, df=df) 
        x <- t(apply(x, 1, function(y) y + gamma_mean[,i]))
        x <- cbind(model$mu$gamma[1, i], x)
        colnames(x) <- rownames(model$mu$gamma)
        return(x)
      })
    }else{
      sim_value <- lapply(1:length(gamma_mean), function(i) {
        invShape <- invV*error.prec[i]
        shape <- solve(invShape)
        x <- rmvt(nsims, sigma=shape, df=df) 
        x <- x + gamma_mean[i]
        x <- cbind(model$mu$gamma[1, i], x)
        colnames(x) <- rownames(model$mu$gamma)
        return(x)
      })
    }

  }
  else if(model$settings$gamma$mode == "LinearRegression"){
    browser()
    gamma_mean <- model$mu$gamma
    rownames(gamma_mean) <- colnames(model$settings$covariates$X)
    gamma_variance <- (model$mu$std.gamma)^2
    sim_value <- lapply(1:ncol(gamma_mean), function(i) {
      x <- mvrnorm(n = nsims, mu = gamma_mean[, i], Sigma = diag(gamma_variance[,i]))
      colnames(x) <- rownames(gamma_mean)
      return(x)
    })
  }
  return(sim_value)
}

calc_relative_change <- function(model, gamma_sim, category = NULL){

  K <- model$settings$dim$K
  # cov <- colnames(scSTMobj$settings$covariates$X)[-1]
  X <- unique(scSTMobj$settings$covariates$X)
  rownames(X) <- apply(X, 1, function(row) {
    case_when(
      row["TimeTime2"] == 1 & row[category] == 1 ~ "alt_t2",
      row["TimeTime2"] == 1 ~ "ref_t2",
      row[category] == 1 ~ "alt_t1",
      TRUE ~ "ref_t1")
  })

  # find the category, if there is any
  category_verify <- colnames(X)[!grepl("(:|Time2)", colnames(X))]

  if(!is.null(category)){
    if (!category %in% category_verify) stop("category needs to be one of ", category_verify)
  }

  results <- list()
  for (i in 1:length(gamma_sim)) {
    results[[i]] <- X %*% t(gamma_sim[[i]])
  }
  
  n_rows <- nrow(results[[1]])
  n_cols <- ncol(results[[1]])
# re-organize the previous output
  result <- vector("list", n_rows)
  for (i in 1:n_cols) {
    combined_columns <- do.call(cbind, lapply(results, function(mat) mat[,i]))
    result[[i]] <- combined_columns
  }
  

  map_to_simplx <- function(x) {
    exp(x - log(sum(exp(x))))
  }
  
  result <- lapply(result, function(mat) cbind(mat, 0))
  simplx_results <- lapply(result, function(mat) {
    x <- t(apply(mat, 1, map_to_simplx)) 
    colnames(x) <- paste0("topic_", 1:K)
    return(x)
  })
  

  # else if(model$settings$gamma$mode == "LinearRegression"){
  #   
  #   gamma_mean <- model$mu$gamma
  #   rownames(gamma_mean) <- sub("^X", "", rownames(gamma_mean))
  #   gamma_variance <- (model$mu$std.gamma)^2
  #   sim_value <- lapply(1:ncol(gamma_mean), function(i) {
  #     x <- mvrnorm(n = nsims, mu = gamma_mean[, i], Sigma = diag(gamma_variance[,i]))
  #     colnames(x) <- rownames(gamma_mean)
  #     return(x)
  #   })
  #   
  #   K <- model$settings$dim$K
  #   cov <- colnames(scSTMobj$settings$covariates$X)
  #   X <- scSTMobj$settings$covariates$X
  #   
  #   # find the category, if there is any
  #   category_verify <- cov[!grepl("(:|Time2)", cov)]
  #   
  #   if(!is.null(category)){
  #     if (!category %in% category_verify) stop("category needs to be one of ", category_verify)
  #   }
  #   
  #   results <- list()
  #   for (i in 1:length(sim_value)) {
  #     if (!is.null(category)) {
  #       # reference group for t1
  #       t1_ref <- sim_value[[i]][,"(Intercept)"]
  #       # reference group for T2
  #       t2_ref <- apply(sim_value[[i]], 1, function(x) sum(x[c("(Intercept)", "TimeTime2")]))
  #       
  #       # t1 for alt group
  #       t1_alt <- apply(sim_value[[i]], 1, function(x) sum(x[c("(Intercept)", category)]))
  #       # t2 for alt group
  #       t2_alt <- apply(sim_value[[i]], 1, function(x) sum(x[c("(Intercept)", "TimeTime2", category, paste0("TimeTime2:", category))]))
  #       # t1_ref <- rep(0, nsims)
  #       # # reference group for T2
  #       # t2_ref <- apply(sim_value[[i]], 1, function(x) sum(x[ "TimeTime2"]))
  #       # 
  #       # # t1 for alt group
  #       # t1_alt <- apply(sim_value[[i]], 1, function(x) sum(x[category])) 
  #       # # t2 for alt group
  #       # t2_alt <- apply(sim_value[[i]], 1, function(x) sum(x[c("TimeTime2", category, paste0("TimeTime2:", category))]))
  #       # Store the results for this category and topic K
  #       results[[paste0("K", i)]] <- list(t1_ref = t1_ref, t2_ref = t2_ref, t1_alt = t1_alt, t2_alt = t2_alt)
  #       sublist_names <- c("t1_ref", "t2_ref", "t1_alt", "t2_alt")
  #     } else {
  #       # If no category, calculate for Time2 vs Time1 only for topic K
  #       t1 <- sim_value[[i]][,"(Intercept)"]
  #       t2 <- apply(sim_value[[i]], 1, function(x) sum(x[c("(Intercept)", "TimeTime2")]))
  #       results[[paste0("K", i)]] <- list(t1 = t1, t2 = t2)
  #       sublist_names <- c("t1", "t2")
  #     }
  #   }
  #   
  #   map_to_simplx <- function(x) {
  #     exp(x - log(sum(exp(x))))
  #   }
  #   
  #   combined_results <- list()
  #   for (name in sublist_names) {
  #     # Use lapply to extract and cbind the corresponding matrices across K1 to K4
  #     combined_matrix <- do.call(cbind, lapply(results, `[[`, name))
  #     combined_matrix <- cbind(combined_matrix, 0)
  #     colnames(combined_matrix) <- paste0("K", 1:ncol(combined_matrix))
  #     combined_results[[name]] <- combined_matrix
  #   }
  #   
  #   simplx_results <- lapply(combined_results, function(mat) {
  #     t(apply(mat, 1, map_to_simplx)) 
  #   })
  # }
  # 
  return(simplx_results)
}
