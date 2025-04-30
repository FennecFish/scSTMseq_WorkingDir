# This script is to run propeller with seurat output
# with sample block in limma
setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
library(speckle)
library(limma)
library(CellBench)
library(mclust)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
seurat_dir <- args[2]
file_name <- args[3]
seurat_name <- file_name

file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat("Seurat dir is",seurat_dir, "\n")
cat(file_name, "\n")

set_level <- sub("^seurat_(.*)\\.rds$", "\\1", file_name)
sims_name <-  paste0(dir, "/sims/sims_",set_level,".rds")

get_transformed_data <- function(seurat.sims, transform="logit"){
  clusters <- Idents(seurat.sims)
  sample <- seurat.sims$Sample
  group <- seurat.sims$Time
  tab <- table(sample, clusters, group)
  props <- apply(tab, c(1, 3), function(x) x / sum(x))
  props <- aperm(props, c(2, 1, 3)) 
  
  if(transform=="asin"){
    message("Performing arcsin square root transformation of proportions")
    prop.trans <- asin(sqrt(props))
  }else if(transform=="logit"){
    message("Performing logit transformation of proportions")
    props.pseudo <- apply(tab, c(1, 3), function(x) (x + 0.5) / sum(x + 0.5))
    props.pseudo <- aperm(props.pseudo, c(2, 1, 3)) 
    prop.trans <- log(props.pseudo/(1-props.pseudo))
  }
  prop.trans <- rbind(prop.trans[, , 1], prop.trans[, , 2])
  colnames(prop.trans) <- paste0("Cluster", 1:ncol(prop.trans))
  prop.trans <- t(prop.trans)
  return(prop.trans)
}

update_propeller <- function(prop.trans, trend=FALSE, robust=TRUE){
  sample <- colnames(prop.trans)
  group <- rep(c(0, 1), each = length(unique(sample)))
  design<-model.matrix(~0+factor(group))
  colnames(design) <- c("Time1", "Time2")
  corfit<-duplicateCorrelation(prop.trans,design,block=sample)
  fit<-lmFit(prop.trans,design, block=sample, correlation=corfit$consensus)
  contrasts <- makeContrasts(Time = Time2 - Time1, levels = design)
  fit.cont <- contrasts.fit(fit, contrasts = contrasts) 
  fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)
  out <- data.frame(coefficients = fit.cont$coefficients[,1],
                    t = fit.cont$t[,1],
                    df.total = fit.cont$df.total,
                    p.value = fit.cont$p.value[,1])
  rownames(out) <- rownames(fit.cont$coefficients)
  out$fdr <- p.adjust(out$p.value, method="BH")
  return(out)
}
################################################################################
############################## read in seurat ##################################
################################################################################
# select ones with the best clustering results
sims <- readRDS(sims_name)
seurat.sims <- readRDS(seurat_name)
smeta <- seurat.sims@meta.data %>% as.data.frame()
sub_sims <- sims[,rownames(smeta)] # filter by the rows
seurat.adj <- sapply(smeta[,4:7], function(x) {
  adjustedRandIndex(x, sub_sims$Group)
})
# select the resolution that has the highest ARI. When there are multiple, select the first one
best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
ARI <- c(basename(dir), set_level, max(seurat.adj))
saveRDS(ARI, paste0("res/1MultiSample_SingleResponse_Simulation/seurat_ARI/ARI_",basename(dir), "_", set_level, ".rds"))

seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
Idents(seurat.sims) = seurat_cluster[,1]
seurat.sims$Time <- colData(sims)$Time[match(rownames(seurat.sims@meta.data),colData(sims)$Cell)]
seurat.sims$Sample <- colData(sims)$Sample[match(rownames(seurat.sims@meta.data),colData(sims)$Cell)]

prop.trans.logit <- get_transformed_data(seurat.sims = seurat.sims, transform="logit")
out.logit <- update_propeller(prop.trans.logit)

dir_path = paste0(dir, "/propeller_limma_sample_block/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "logit_", set_level, ".rds")
saveRDS(out.logit, save_file_name)

prop.trans.asin <- get_transformed_data(seurat.sims = seurat.sims, transform="asin")
out.asin <- update_propeller(prop.trans.asin)

save_file_name <- paste0(dir_path, "asin_", set_level, ".rds")
saveRDS(out.asin, save_file_name)
