#' Subset the stm results based on sample ID or a subset of document names
#' 
#' @param stmobj An \code{STM} object created by \code{\link{multi_stm}}
#' @param subDocName A vector of document names in the STM object 
#' @param sampleID A vector of sample IDs that used in the \code{\link{multi_stm}}
#' @export
STMsubset <- function(stmobj, subDocName = NULL, sampleID = NULL){
    
    if(is.null(subDocName) & is.null(sampleID)) stop("At least one of document sets or sample ID of interet need to be specified")
    if(!is.null(sampleID) & is.null(subDocName)) {
        subDocName <- stmobj$DocName[stmobj$sampleID %in% sampleIDs]
    }
    if(!is.null(sampleID) & !is.null(subDocName)) {
        subDocName <- union(stmobj$DocName[stmobj$sampleID %in% sampleIDs], subDocName)
    }
    
    # booled doc index
    doc.index <- stmobj$DocName %in% subDocName
    stmobj$DocName <- stmobj$DocName[doc.index]
    stmobj$sampleID <- stmobj$sampleID[doc.index]
    # subset mu
    stmobj$mu$mu <- stmobj$mu$mu[,doc.index]
    
    #change settings
    stmobj$settings$dim$N <- sum(doc.index)
    stmobj$settings$dim$I <- length(unique(stmobj$sampleID))
    stmobj$settings$dim$wcounts <- NULL
    stmobj$settings$dim$samples <- stmobj$sampleID
    stmobj$settings$covariates$X <- stmobj$settings$covariates$X[doc.index,]
    stmobj$settings$covariates$betaindex <- stmobj$settings$covariates$betaindex[doc.index]
    
    # other variables
    stmobj$theta <- stmobj$theta[doc.index,]
    stmobj$eta <- stmobj$eta[doc.index,]
    return(stmobj)
}