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

# sims <- scSTMobj$settings$sce
# dir_path <- paste0(dir, "/1000sims/sims_", set_level, ".rds")
# saveRDS(sims, dir_path)
Manova.res <- MaxThetaManova(model = scSTMobj, nsims = 100, Sample = "Sample", Time = "Time")

thetaManova_dir <- sub("scSTM", "ManovaTheta", scSTM_dir)
dir_path <- paste0(dir, "/", thetaManova_dir)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
saveRDS(Manova.res, file = paste0(dir_path, "/ManovaRM_Max_",set_level,".rds"))

# rm.wts <- lapply(Manova.res, function(x){
#   # use the Manova with repeated measure for power and Type I error
#   x <- x$Interaction
#   if(x$Warnings == "No Warnings"){
#     dat <- x$resampling[1]
#   }else{
#     dat = NA
#   }
#   return(dat)
# })

rm.mats <- lapply(Manova.res, function(x){
  dat <- x$resampling[2]
  return(dat)
})

dat.clean <- function(given_list){
  given_list <- do.call(rbind, given_list)
  rownames(given_list) <- paste0("replicate", 1:nrow(given_list))
  colnames(given_list) <- "pValue"
  return(given_list)
}
res <- list(rm.mats = dat.clean(rm.mats), sce = scSTMobj$settings$sce)
saveRDS(res, file = paste0(dir, "/", thetaManova_dir, "/Manova_Max_pValue_",set_level,".rds"))
