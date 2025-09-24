setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)
library(scater)
library(scran)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
file_name <- args[2]
design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
cat("Design is", design, "\n")
file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat(file_name, "\n")

# set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
# sims <- readRDS(paste0(dir, "/sims/", file_name))
set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
file_name <- paste0("sims_", set_level, ".rds")
sims <- readRDS(paste0(dir, "/sims/", file_name))

seurat.sims <- CreateSeuratObject(counts = counts(sims))
# seurat.sims <- subset(seurat.sims, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seurat.sims <- NormalizeData(seurat.sims)
# Visualize QC metrics as a violin plot
# VlnPlot(seurat.sims, features = c("nFeature_originalexp", "nCount_originalexp"), ncol = 2)
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(seurat.sims)
seurat.sims <- ScaleData(seurat.sims, features = all.genes)
seurat.sims <- RunPCA(seurat.sims, features = VariableFeatures(object = seurat.sims))

seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
seurat.sims <- FindClusters(seurat.sims, resolution = seq(0.5,2,0.5))

dir_path <- paste0(dir, "/seurat/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "seurat_", set_level, ".rds")
saveRDS(seurat.sims, file = save_file_name)