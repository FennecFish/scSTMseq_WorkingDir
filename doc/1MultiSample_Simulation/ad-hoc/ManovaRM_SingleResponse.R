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

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
scSTM_dir <- args[2]
file_name <- args[3]
file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat("scSTM dir is",scSTM_dir, "\n")
cat(file_name, "\n")

set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", file_name)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

scSTM_name <- paste0(dir, "/", scSTM_dir, "/scSTM_",
                     set_level,".rds")
if(file.exists(scSTM_name)){
  scSTMobj <- readRDS(scSTM_name)
  scSTMobj <- select_top_scSTM(scSTMobj)
}else{
  next
}

PosteriorRep <- PosteriorPropRep(model = scSTMobj, nsims = 100, Sample = "Sample", Time = "Time", Group = NULL)
Manova.res <- ThetaManova(PosteriorRep, Time = "Time", Group = NULL, use_mean = F)

thetaManova_dir <- "ManovaRM"
dir_path <- paste0(dir, "/", thetaManova_dir)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
saveRDS(ThetaManova, file = paste0(dir_path, "/ManovaRM_",set_level,".rds"))

rm.wts <- lapply(Manova.res, function(x){
  # use the Manova with repeated measure for power and Type I error
  x <- x$Interaction
  if(x$Warnings == "No Warnings"){
    dat <- x$resampling[1]
  }else{
    dat = NA
  }
  return(dat)
})

rm.mats <- lapply(Manova.res, function(x){
  x <- x$Interaction
  dat <- x$resampling[2]
  return(dat)
})

dat.clean <- function(given_list){
  given_list <- do.call(rbind, given_list)
  rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
  colnames(given_list) <- "pValue"
  return(given_list)
}
res <- list(rm.wts = dat.clean(rm.wts), rm.mats = dat.clean(rm.mats), sce = scSTMobj$settings$sce)
saveRDS(res, file = paste0(dir, "/", thetaManova_dir, "/Manova_pValue_",set_level,".rds"))

# sims <- scSTMobj$settings$sce
# # collapse simulate cells by their group and sample
# theta.collapsed <- colData(sims) %>%
#   as.data.frame() %>%
#   group_by(Time, Sample, Group, Response) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   group_by(Time, Sample) %>%
#   dplyr::mutate(proportion = count / sum(count)) %>%
#   ungroup() %>%
#   dplyr::select(-count) %>%
#   pivot_wider(names_from = Group, values_from = proportion, values_fill = 0)
# 
# meta <- theta.collapsed %>% dplyr::select(!starts_with("Group"))
# cluster <- theta.collapsed %>% dplyr::select(starts_with("Group"))
# 
# cluster <- compositions::ilr(cluster)
# colnames(cluster) <- paste0("Group", 1:ncol(cluster))
# data <- cbind(meta, cluster)
# ##### Fit Manova.RM #######
# response_vars <- grep("^Group", colnames(data), value = TRUE)
# fit <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Time")),
#               data = data, subject = "Sample", within = "Time", iter = 1000)
# saveRDS(fit, file = paste0(dir, "/", thetaManova_dir, "/Manova_sims_",set_level,".rds"))
# cat("Complete SIM Manova \n")
