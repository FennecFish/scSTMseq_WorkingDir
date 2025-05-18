setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(stringr)
library(ggplot2)
library(dplyr)

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test4/", pattern = "scSTM*")
# res <- vector(mode="list",length=length(files))
res <- data.frame()
for (file in files){
  set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file)
  alt <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test4/", 
                        file))
  if(class(alt) == "selectModel"){
    all_values <- unlist(alt$bound)
    max_value <- max(all_values)
    max_position_in_vector <- which(all_values == max_value)
    alt <- alt$runout[[max_position_in_vector]]
  }
  K <- alt$settings$dim$K
  file.res <- data.frame()
  # sample from the normal
  for (k in 1:(K-1)){
    X <- rnorm(100, mean = alt$mu$gamma[2,k], sd = alt$mu$sn[k])
    tres <- t.test(X, mu = 0)
    temp <- c(unlist(str_split(set_level, "_")), k, tres$statistic, tres$p.value)
    file.res <- rbind(file.res, temp)
  }
  colnames(file.res) <- c("seed", "control", "level", "nCellType", "k", "statistic", "p.value")
  res <- rbind(res, file.res)
}
write.csv(res, file = "res/cc_single_ttest.csv")

res <- read.csv("res/cc_single_ttest.csv")
