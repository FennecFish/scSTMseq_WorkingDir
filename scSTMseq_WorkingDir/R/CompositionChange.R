composition_change <- function(stmobj){
    time <- stmobj$settings$covariates$X[,-1]
    dat <- stmobj$theta
    rownames(dat) <- stmobj$DocName
    t1 <- dat[which(time==1),]
    t2 <- dat[which(time==2),]
    names(stmobj$sampleID) <- stmobj$DocName
    K <- stmobj$settings$dim$K
    sampleID <- unique(stmobj$sampleID)
    res <- data.frame()
    for(sample in sampleID){
        x1 <- t1[rownames(t1) %in% names(which(stmobj$sampleID==sample)),]
        x2 <- t2[rownames(t2) %in% names(which(stmobj$sampleID==sample)),]
        
        y1 <- compositions::ilr(x1)
        y2 <- compositions::ilr(x2)
        m1 <- colMeans(y1)
        m2 <- colMeans(y2)
        d <- dim(y1)[2]
        n1 <- dim(y1)[1] 
        n2 <- dim(y2)[1]
        s1 <- (crossprod(y1) - n1 * tcrossprod(m1) ) / n1
        s2 <- (crossprod(y2) - n2 * tcrossprod(m2) ) / n2
        
        i <- 1
        s1h <- s1
        s2h <- s2
        s1inv <- solve(s1h)  ;  s2inv <- solve(s2h)
        mha <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2)
        s1h <- s1h + tcrossprod(m1 - mha)
        s2h <- s2h + tcrossprod(m2 - mha)
        s1inv <- solve(s1h) 
        s2inv <- solve(s2h)
        mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
        while ( sum( abs(mha - mhb) ) > 1e-6 ) {
            i <- i + 1
            mha <- mhb
            mhb <- solve( n1 * s1inv + n2 * s2inv, n1 * s1inv %*% m1 + n2 * s2inv %*% m2 )
            s2h <- s1h + tcrossprod(m1 - mhb)
            s2h <- s1h + tcrossprod(m2 - mhb)
        }
        dof <- d
        stat <- n1 * log( det(s1h) / det(s1) ) + n2 * log( det(s2h) / det(s2) )
        pvalue <- pchisq(stat, dof, lower.tail = FALSE)
        
        # summarize all data
        res <- rbind(res, c(stat,pvalue,dof))
    }
    rownames(res) <- sampleID
    colnames(res) <- c("Statistics", "p-value","DegreeOfFreedom")
    # res <- res[,c(1,2,4)]
    return(res)
}

plot_dist <- function(stmobj, sample){
    time <- stmobj$settings$covariates$X[,-1]
    dat <- stmobj$theta
    rownames(dat) <- stmobj$DocName
    t1 <- dat[which(time==1),]
    t2 <- dat[which(time==2),]
    names(stmobj$sampleID) <- stmobj$DocName
    x1 <- t1[rownames(t1) %in% names(which(stmobj$sampleID==sample)),]
    x2 <- t2[rownames(t2) %in% names(which(stmobj$sampleID==sample)),]
    y1 <- compositions::ilr(x1)
    y2 <- compositions::ilr(x2)
    d1 <- dist(y1)
    d2 <- dist(y2)
    
    return(boxplot(d1,d2))
}

process_scSTM <- function(scSTMobj, theta = NULL) {
  if(is.null(theta)){theta = scSTMobj$theta} else{theta = theta}
  max_indices <- apply(theta, 1, which.max)
  colnames(theta) <- paste0("topic_", 1:ncol(theta))
  rownames(theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(theta)[max_indices]
  names(res_cluster) <- rownames(theta)
  return(res_cluster)
}

propeller_change <- function(stmobj, sample, nsims){
  cluster <- vector(mode = "list", length = nsims)
  thetasims <- thetaPosterior(stmobj, nsims=nsims, type="Global")
  new_lists <- vector("list", nsims)
  for (i in 1:nsims) {
    new_lists[[i]] <- lapply(thetasims, function(x) x[i, ])
  }
  new_matrices <- lapply(new_lists, function(lst) {
    do.call(rbind, lst)
  })
  
  thetasims <- new_matrices
  cluster <- lapply(thetasims, function(mat) {
    process_scSTM(stmobj, theta = mat)
  })

  vec.cluster <- unlist(cluster)
  grp <- rep(stmobj$settings$covariates$X[,2], nsims)
  samp <- rep(1:nsims, each = length(stmobj$DocName))
  p <- propeller(clusters = vec.cluster, group = grp, sample = samp, trend = FALSE, robust = TRUE)
  return(p)
}


