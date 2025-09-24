setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(scater)
library(scran)
library(Seurat)

dat <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/single_patient/BIOKEY_10.rds")

dat <- logNormCounts(dat)
logcounts(dat) <- as.matrix(logcounts(dat))
dat <- runPCA(dat)
plotPCA(dat, colour_by = "cellType")
plotPCA(dat, colour_by = "timepoint")
# seurat <- as.Seurat(dat, counts = "counts", data = "logcounts")
rm(dat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
DimPlot(seurat, reduction = "pca")
        
seurat <- RunUMAP(seurat, dims = 1:10)
saveRDS(seurat, "/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/anti_PD1_cohort1_seurat.rds")

DimPlot(pbmc, reduction = "umap")

logcounts(dat) <- as.matrix(logcounts(dat))

