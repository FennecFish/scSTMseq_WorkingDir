convert_count_matrix <- function(count_matrix) {
    
    gene_indices <- seq_len(ncol(count_matrix))
    
    cell_list <- apply(count_matrix, 1, function(cell_counts) {
        # Identify nonzero gene counts
        nonzero_indices <- which(cell_counts != 0)
        
        # Filter out zeros from counts and corresponding gene indices
        if (length(nonzero_indices) > 0) {
            nonzero_gene_indices <- gene_indices[nonzero_indices]
            nonzero_counts <- cell_counts[nonzero_indices]
            
            # Construct the matrix with filtered indices and counts
            mat <- rbind(nonzero_gene_indices, nonzero_counts)
            dimnames(mat) <- NULL
            return(list(mat))
        } 
    })
    
    # Simplify the list structure if nested due to list conversion in apply
    cell_list <- unlist(cell_list, recursive = FALSE)
    
    # Set names for each list element
    names(cell_list) <- rownames(count_matrix)
    return(cell_list)
}


# require SingleCellExperiment, with meta data in ColData
prepsce <- function(sce, sample = NULL, lower.thresh=1, upper.thresh=Inf, 
                    verbose=TRUE, filter = TRUE){
    #Functions:
    # 1) Detect and renumber zero-indexed data.
    # 2) Detect and renumber missing terms
    # 3) Remove words appearing only in [lower.thresh] documents.
    # 4) Remove words appear in upper.thresh or more of the documents.
    
    # Can also optionally subsample the data.
    
    # check if used SingleCellExperiment
    if((!class(sce)[1]=="SingleCellExperiment")){
        stop("The input data needs to be a SingleCellExperiment Object")
    }

    gene_sums <- colSums(counts(sce))
    sce <- sce[, gene_sums > 0]
    if(sum(gene_sums==0) > 0) {
        cat("Removing", sum(gene_sums==0), "genes have 0 cells...\n")
    } 
    
    cell_sums <- rowSums(counts(sce))
    sce <- sce[cell_sums > 0, ]
    if(sum(cell_sums==0) > 0) {
        cat("Removing", sum(cell_sums==0), "cells have 0 genes...\n")
    } 
    
    # extract documents (cell) and vocab (genes) 
    documents <- convert_count_matrix(t(assays(sce)$counts))

    vocab <- rownames(sce)
    
    #error check for inputs
    if((is.null(documents))){
        stop("No cells found in your data")
    }
    if(is.null(vocab)){
        stop("No genes found in your data")
    }
    if(!inherits(documents,"list")) {
        stop("documents must be a list in stm() format.  See ?stm() for format.  
          See ?readCorpus() for tools for converting from popular formats")
    }
    # 
    # if(!is.null(subsample)) {
    #     index <- sample(1:length(documents), subsample)
    #     documents <- documents[index]
    #     if(!is.null(meta)) meta <- meta[index, , drop = FALSE] 
    # }
    # browser()
    #check that there are no 0 length documents

    
    # triplet <- doc.to.ijv(documents) #this also fixes the zero indexing.
    # nms <- names(documents) 
    # documents <- ijv.to.doc(triplet$i, triplet$j, triplet$v)
    # names(documents) <- nms
    meta <- colData(sce)
    docs.removed <- c()
    
    #Detect Missing Terms
    miss.vocab <- NULL
    vocablist <- sort(seq_along(rownames(sce)))
    nonzero_counts <- t(counts(sce)) > 0
    wordcounts <- colSums(nonzero_counts)
    names(wordcounts) = NULL

    if(length(vocablist)<length(vocab)) {
        if(verbose) cat("Detected Missing Terms, renumbering \n")
        miss.vocab <- vocab[-vocablist]
        vocab <- vocab[vocablist]
        new.map <- cbind(vocablist, 1:length(vocablist))
        documents <- lapply(documents, function(d) {
            nm <- names(d)
            d[1,] <- new.map[match(d[1,], new.map[,1]),2]
            names(d) <- nm
            return(d)
        })
        wordcounts <- wordcounts[vocablist]
    }
    
    #Remove Words Appearing Only n Times
    if(filter){
      toremove <- which(wordcounts <= lower.thresh | wordcounts >= upper.thresh)
      keepers <- which(wordcounts > lower.thresh & wordcounts < upper.thresh)
      droppedwords <- c(miss.vocab,vocab[toremove])
    }else{
      toremove <- as.integer()
      keepers <- seq_along(wordcounts)
      droppedwords <- c(miss.vocab,vocab[toremove])
    }
    
    if(length(toremove)) {
        if(verbose) cat(sprintf("Removing %i of %i genes (%i of %i tokens) due to frequency \n", 
                                length(toremove), length(wordcounts), sum(wordcounts[toremove]), sum(wordcounts)))
        vocab <- vocab[-toremove]
        remap <- 1:length(keepers)
        for(i in 1:length(documents)) {
            doc <- documents[[i]]
            dockeep <- doc[1,]%in%keepers
            doc <- doc[,dockeep,drop=FALSE]
            doc[1,] <- remap[match(doc[1,], keepers)]
            documents[[i]] <- doc
            if(ncol(doc)==0) docs.removed <- c(docs.removed,i)
        }
        if(length(docs.removed)) {
            if(verbose) cat(sprintf("Removing %i Cells with No Gene Counts", length(docs.removed), "\n"))
            documents <- documents[-docs.removed]
        }
        toprint <- sprintf("Your Dataset now has %i cells, %i genes and %i tokens. \n", 
                           length(documents), length(vocab), sum(wordcounts[keepers]))
        if(verbose) cat(toprint)
    }
    
    # re-organize meta data to remove cells that has been removed
    if(!is.null(docs.removed) & !is.null(meta)){
        meta <- meta[-docs.removed, , drop = FALSE]
    }
    # re-organize sce to remove vocab and docs
    sce <- sce[vocab, rownames(meta)]
    #recast everything as an integer
    documents <- lapply(documents, function(x) matrix(as.integer(x), nrow=2))
    if(!is.null(sample)) {
        sample_list <- meta[sample]
        return(list(documents=documents, vocab=vocab, 
                    sample = sample_list, meta=meta, sce=sce,
                    words.removed=droppedwords, docs.removed=docs.removed, 
                    tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
    } else {
        return(list(documents=documents, vocab=vocab, 
                    sample = NULL, meta=meta, sce=sce, 
                    words.removed=droppedwords, docs.removed=docs.removed, 
                    tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
    }
    
    
   
}
