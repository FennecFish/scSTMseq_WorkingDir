setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
library(speckle)
library(limma)
library(mclust)

#################################################################################
################################# ARI #########################################
#################################################################################
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
# design <- "nSample20_nCellType5_noBatch_StromalCell/"
# scSTMdir <- "scSTM_Pooled_Content_Sample_Prevalence_Time"
scSTMdir <- "scSTM_Pooled_noContent_Prevalence_Time"

# extract scSTM files from different design
paths <- Sys.glob(file.path(dir, "*nSample5*_noBatch_CancerCell*", scSTMdir)) 
files <- unlist(lapply(paths, list.files, full.names = TRUE))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

dat <- data.frame()
for (file in ari_file){
  temp <- readRDS(paste0("res/1MultiSample_SingleResponse_Simulation/seurat_ARI/", file))
  dat <- rbind(dat, temp)
}
colnames(dat) <- c("design", "set_level", "ARI")
res.temp <- data.frame(
  modelType = sapply(strsplit(dat$set_level, "_"), `[`, 2),
  seed = sapply(strsplit(dat$set_level, "_"), `[`, 1),
  nCellType = as.numeric(gsub("nCellType", "", sapply(strsplit(dat$design, "_"), `[`, 2))),
  nSample = as.numeric(gsub("nSample", "", sapply(strsplit(dat$design, "_"), `[`, 1))),
  Batch = ifelse(sapply(strsplit(dat$design, "_"), `[`, 3)=="noBatch", FALSE, TRUE),
  CancerType = ifelse(sapply(strsplit(dat$design, "_"), `[`, 4)=="StromalCell", FALSE, TRUE),
  ARI = as.numeric(dat$ARI)
)
write.csv(res.temp, file = "res/1MultiSample_SingleResponse_Simulation/seurat_ARI.csv")


# ari <- read.csv("/proj/milovelab/wu/scLDAseq/res/1MultiSample_SingleResponse_Simulation/seurat_ARI.csv")
# 
# ari <- ari %>%
#   as.data.frame() %>%
#   dplyr::mutate(scSTMTypeNumeric = case_when(
#     modelType == "NullModel" ~ 0,
#     grepl("HighVar", modelType) ~ as.numeric(sub("HighVar", "", modelType))
#   ))
# 
# ggplot(ari, aes(x = factor(scSTMTypeNumeric), y = ARI)) +
#   geom_boxplot() +
#   labs(
#     x = "scSTMType (0: NullModel, SD values)",
#     y = "ARI",
#     title = "ARI for Clustering Accuracy Using Seurat"
#   ) +
#   theme_minimal()
