make.heldout.sce <- function(sce , N=floor(.1*ncol(sce)), 
                         proportion=.5, seed=NULL, sample = NULL) {
  if(!is.null(seed)) set.seed(seed)
  
    if(is.null(sce)) stop("Please provide a SingleCellExperiment Object")
    
    # Convert the corpus to the internal STM format
    # args <- prepsce(sce)
    # documents <- args$documents
    # vocab <- args$vocab
    
    ncell <- ncol(sce)
    index <- sort(sample(1:ncell, N))
    pie <- proportion
    
    # in held-out set, select 50% of genes to take out
    missing <- vector(mode="list", length=N)
    missing.gene.index <- vector(mode="list", length=N)
    # missing$index <- index # this is the cells that are missing
    
    training.sce <- sce
    ct <- 0
    for(i in index){
        ct <- ct + 1
        gene_counts <- counts(training.sce)[,i]
        nsamp <- floor(pie*sum(gene_counts))
        missing.index <- sample(1:sum(gene_counts),nsamp)
        
        gene_vector <- rep(names(gene_counts), times = gene_counts)
        missing_genes <- gene_vector[missing.index]
        missing_gene_counts <- table(missing_genes)
        missing[[ct]] <- missing_gene_counts
        
        train_gene_counts <- gene_counts
        m.gene.index <- match(names(missing_gene_counts),names(gene_counts))
        train_gene_counts[m.gene.index] <-
            gene_counts[m.gene.index] - missing_gene_counts      
        
        counts(training.sce)[,i] <- train_gene_counts
        missing.gene.index[[ct]] <- m.gene.index
        
    }
    
    missing <- list(index=index, cells=missing, missing.gene.index = missing.gene.index)
    heldout <- list(training.sce = training.sce, heldout.dat=missing)
    class(heldout) <- "heldout"
    return(heldout)
  #   # use whichever give more samples N or proportion
  #   ndoc <- ncol(sce)
  #   pie <- ifelse(proportion * ndoc >= N, proportion, N/ndoc)
  #   
  #   # Split the meta data by sample ID
  #   if(!is.na(sample)){
  #       groups <- split(seq_along(sce[[sample]]), sce[[sample]])
  #       
  #       # Sample from each group certain proportion
  #       sampled_indices <- lapply(groups, function(x) 
  #       {sample(x, size = floor(length(x) * pie))})
  #       sampled_indices <- sort(unlist(sampled_indices, use.names = FALSE))
  #   } else{
  #       # if no sample provided
  #       sampled_indices <- sort(sample(seq_along(1:ncol(sce)), size = floor(ncol(sce) * pie)))
  #   }
  #   
  #   sce_sub <- sce[,sampled_indices]
  #   missing <- seq_along(1:ncol(sce))[!seq_along(1:ncol(sce)) %in% sampled_indices]
  # 
  # return(list(sce_sub = sce_sub, keep = sampled_indices, missing = missing))
}

#' #' @export
#' eval.heldout <- function(model, missing) {
#'   heldout <- vector(length=length(missing$index))
#'   ntokens <- vector(length=length(missing$index))
#'   beta <- lapply(model$beta$logbeta, exp)
#'   bindex <- model$settings$covariates$betaindex[missing$index]
#'   for(i in 1:length(missing$index)) {
#'     docid <- missing$index[i]
#'     words <- missing$docs[[i]][1,]
#'     probs <- model$theta[docid,]%*%beta[[bindex[i]]][,words]
#'     probs <- rep(probs, missing$docs[[i]][2,])
#'     heldout[i] <- mean(log(probs)) 
#'     ntokens[i] <- sum(missing$docs[[i]][2,])
#'   }
#'   out <- list(expected.heldout=mean(heldout, na.rm=TRUE), doc.heldout=heldout,
#'               index=missing$index, ntokens=ntokens) #the mean here to deal with 0 length docs
#'   return(out)
#' }