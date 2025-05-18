setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
require(pals)
library(MASS)
library(tibble)
library(MANOVA.RM)
library(mclust)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
suerat_dir <- args[2]
seurat_name <- args[3]
file_name <- basename(seurat_name)
cat("Dir is", dir, "\n")
cat("seurat dir is",suerat_dir, "\n")
cat(file_name, "\n")

set_level <- sub("^seurat_(.*)\\.rds$", "\\1", file_name)
sims_name <-  paste0(dir, "/sims/sims_",set_level,".rds")

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

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
best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]

seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% dplyr::select(all_of(best_res))
Idents(seurat.sims) = seurat_cluster[,1]

dataset <- data.frame(
  Cluster = paste0("Cluster", Idents(seurat.sims)),
  Time = colData(sims)$Time[match(rownames(seurat.sims@meta.data),colData(sims)$Cell)],
  Sample = colData(sims)$Sample[match(rownames(seurat.sims@meta.data),colData(sims)$Cell)]
)

theta.collapsed <- dataset %>%
  as.data.frame() %>%
  group_by(Time, Sample, Cluster) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Time, Sample) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  dplyr::select(-count) %>%
  pivot_wider(names_from = Cluster, values_from = proportion, values_fill = 0)

meta <- theta.collapsed %>% dplyr::select(!starts_with("Cluster"))
cluster <- theta.collapsed %>% dplyr::select(starts_with("Cluster"))

cluster <- compositions::ilr(cluster)
colnames(cluster) <- paste0("Cluster", 1:ncol(cluster))
data <- cbind(meta, cluster)
##### Fit Manova.RM #######
response_vars <- grep("^Cluster", colnames(data), value = TRUE)
fit <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
              data = data, subject = "Sample", within = "Time", iter = 1000)

dir_path <- paste0(dir, "/Manova_Seurat")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
saveRDS(fit, file = paste0(dir, "/Manova_Seurat/Manova_seurat_",set_level,".rds"))

rm.mats <- fit$resampling[,2]
saveRDS(rm.mats, file = paste0(dir, "/Manova_Seurat/Manova_seurat_pValue_",set_level,".rds"))

