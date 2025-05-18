selectModel_parallel <- function(sce , K, sample = NULL,
                        prevalence=NULL, content=NULL,
                        max.em.its=100, verbose=TRUE, # init.type = "TopicScore",
                        emtol= 1e-06, seed=NULL, 
                        ts_runs = 10, random_run = 20, frexw=.7, 
                        net.max.em.its=5, netverbose=TRUE, M=10, N=NULL,
                        to.disk=F, control=list(), gc = 5, ...){

    cl <- makeCluster(gc,outfile="")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
    # Ensure the workers source the necessary files
    clusterEvalQ(cl, {
        library(Rcpp)
        library(slam)
        library(SingleCellExperiment)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(stats)
        library(MASS)
        library(mclust)
        r.file <- paste0("R/", list.files("R/"))
        sapply(r.file, source)
        sourceCpp("src/STMCfuns.cpp")
    })
    
    if(!is.null(seed)) set.seed(seed)
    
    runs <- ts_runs + random_run + 1
    
    if(is.null(N)){
        N <-  round(.2*runs)
    }
    
    if(runs < 2){
        stop("Number of runs must be two or greater.")
    }
    
    if(runs < N){
        stop("Number in the net must be greater or equal to the number of final models.")
    }
    
    args <- prepsce(sce)
    documents <- args$documents
    vocab <- args$vocab
    data <- args$meta
    sce <- args$sce
    # divide runs between using Poisson NMF (Random), TopicScore and Spectral
    # if (runs >=3 ) {
    #     random_run <- ceiling(runs/2)-1
    #     ts_run <- runs-random_run-1
    # }else {
    #     ts_run <- runs - 1
    #     random_run <- 0
    # }
    # 
    
    seedout <- vector("numeric", runs)
    likelihood <- vector("numeric", runs)
    mode <- vector("character", runs)
    
    cat("Casting net \n")

    # Run the models in parallel
    results <- foreach(i = 1:runs, .packages = c('Rcpp'), .combine = 'rbind', .multicombine = TRUE) %dopar% {
      tryCatch({
            if (i <= ts_runs) {
                init_type <- "TopicScore"
            } else if (i == ts_runs + 1) {
                init_type <- "Spectral"
            } else {
                init_type <- "Random"
            }
          message("Running iteration: ", i, " with init_type: ", init_type, "\n")
            mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data,
                                sample = sample, K = K,
                                prevalence = prevalence, content = content, init.type = init_type,
                                max.em.its = net.max.em.its, emtol = emtol, verbose = netverbose, ...)
            message("Completed iteration: ", i, "\n")
            c(seed = mod.out$settings$seed,
              likelihood = mod.out$convergence$bound[length(mod.out$convergence$bound)],
              init_mode = mod.out$settings$init$mode)
        }, error = function(e) {
          message("An error occurred in initiation with", init_type, ": ", e$message, "\n")
            c(seed = NA, likelihood = NA, init_mode = NA)
        })

    }

    results <- results %>% 
        as.data.frame() %>%
        filter(!is.na(likelihood)) %>%
        arrange(likelihood, decreasing = T)
    
    if(dim(results)[1] < N){
        keep <- results
        cat("Return", dim(keep)[1], "models instead. For details, please see error message returned \n")
        N <- dim(results)[1]
    } else{
        keep <- results[1:N,]
    }

    cat("Running select models \n")
    
    runout <- vector("list", N)
    semcoh <- vector("list", N)
    exclusivity <- vector("list", N)
    sparsity <- vector("list", N)
    bound <- vector("list", N)

    final_results <- foreach(i = 1:N, .packages = c('Rcpp'), .combine = 'list', .multicombine = TRUE) %dopar% {
        tryCatch({
            initseed <- as.numeric(keep$seed[i])
            init_type <- keep$init_mode[i]
            #list(initseed= initseed, init_type = init_type)
            message("Running final iteration: ", i, " with init_type: ", init_type, "\n")
            mod.out <- scSTMseq(sce = sce, documents = documents, vocab = vocab, data = data,
                                sample = sample, K = K,
                                prevalence = prevalence, content = content, init.type = init_type,
                                seed = initseed, max.em.its = max.em.its, emtol = emtol,
                                verbose = verbose, ...)
            
            message("Completed final iteration: ", i, "\n")
            # list(runout = mod.out, bound = max(mod.out$convergence$bound))
            # 
            semcoh <- semanticCoherence(mod.out, documents, M)
            if (length(mod.out$beta$logbeta) < 2) {
                exclusivity <- exclusivity(mod.out, M = M, frexw = .7)
                sparsity <- "Sparsity not calculated for models without content covariates"
            } else {
                exclusivity <- "Exclusivity not calculated for models with content covariates"
                kappas <- t(matrix(unlist(mod.out$beta$kappa$params), ncol = length(mod.out$beta$kappa$params)))
                topics <- mod.out$settings$dim$K
                numsparse <- apply(kappas[(K + 1):nrow(kappas), ], 1, function(x) sum(x < emtol))
                sparsity <- numsparse / ncol(kappas)
            }
            
            list(runout = mod.out, bound = max(mod.out$convergence$bound),
                 semcoh = semcoh, exclusivity = exclusivity, sparsity = sparsity)
        }, error = function(e) {
          message("An error occurred in final model with seed", initseed, "and initial type", 
                init_type, ": ", e$message, "\n")
            list(runout = NA, bound = NA, semcoh = NA, exclusivity = NA, sparsity = NA)
        })
    }

    if(N == 1){
        runout <- final_results$runout
        bound <- final_results$bound
        semcoh <- final_results$semcoh
        exclusivity <- final_results$exclusivity
        sparsity <- final_results$sparsity
    } else{
        for (i in 1:N) {
            runout[[i]] <- final_results[[i]]$runout
            bound[[i]] <- final_results[[i]]$bound
            semcoh[[i]] <- final_results[[i]]$semcoh
            exclusivity[[i]] <- final_results[[i]]$exclusivity
            sparsity[[i]] <- final_results[[i]]$sparsity
        } 
    }
    out <- list(runout=runout, bound = bound, semcoh=semcoh, exclusivity=exclusivity, sparsity=sparsity)
    class(out) <- "selectModel"
    return(out)
}