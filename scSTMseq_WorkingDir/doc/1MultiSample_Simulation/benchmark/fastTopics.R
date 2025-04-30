# Goal: Apply fastTopics to datasets for sensitivity analysis 
#       The objective is to test if scSTM has better performance than its intialization point
# Author: Euphy Wu

setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(fastTopics)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
file_name <- args[2]
design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
cat("Design is", design, "\n")
file_name <- basename(file_name)
cat("Dir is", dir, "\n")
cat(file_name, "\n")
gc <- as.integer(args[3])

set_level <- sub("^sims_(.*)\\.rds$", "\\1", file_name)
sims <- readRDS(paste0(dir, "/sims/", file_name))

sims <- sims[rowSums(counts(sims)) != 0,]
counts <- t(counts(sims))
ngroup <- length(unique(sims$Group))
fit <- fit_topic_model(counts,k = ngroup, control.main = list(nc = gc), verbose = "progressbar")

dir_path <- paste0(dir, "/fastTopics/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "fastTopics_", set_level, ".rds")
saveRDS(fit, file = save_file_name)


