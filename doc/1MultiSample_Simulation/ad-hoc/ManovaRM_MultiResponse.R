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
library(compositions)

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


PosteriorRep <- PosteriorPropRep(model = scSTMobj, nsims = 100, Sample = "Sample", Time = "Time", Group = "Response")
Manova.res <- ThetaManova(PosteriorRep, Time = "Time", Group = "Response")

thetaManova_dir <- sub("scSTM", "ManovaTheta", scSTM_dir)
dir_path <- paste0(dir, "/", thetaManova_dir)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
save_file_name <- paste0(dir_path, "/ManovaRM_",set_level,".rds")

thetaManova_dir <- sub("scSTM", "ManovaTheta", scSTM_dir)
saveRDS(Manova.res, file = save_file_name)
cat("Complete Theta Manova \n")

rm.mats <- lapply(Manova.res, function(x){
  # interaction
  interact <- x[["Interaction"]]
  int.dat <- t(interact$resampling[,2]) %>% as.data.frame()
  dat <- int.dat
  dat$Type <- "Interaction"
  # main <- x[["MainEffect"]]
  # if(sum(is.na(main)) == 0){
  #   main.dat <- t(main$resampling[,2]) %>% as.data.frame()
  #   dat <- int.dat %>% 
  #     full_join(main.dat, by = join_by(Response, Time))
  #   dat$Type <- c("Interaction", "Main")
  # }else{
  #   dat <- int.dat
  #   dat$Type <- "Interaction"
  # }
  return(dat)
})

rm.mats <- do.call(rbind, rm.mats)
# rm.mats$Replicates <- rep(paste0("replicate", 1:(nrow(rm.mats)/2)), each =2)
rm.mats$Replicates <- rep(paste0("replicate", 1:nrow(rm.mats)))
saveRDS(rm.mats, file = paste0(dir, "/", thetaManova_dir, "/Manova_pValue_",set_level,".rds"))
cat("Complete Theta Manova MATS \n")

sims <- scSTMobj$settings$sce
# collapse simulate cells by their group and sample
theta.collapsed <- colData(sims) %>%
  as.data.frame() %>%
  group_by(Time, Sample, Group, Response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Time, Sample) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  dplyr::select(-count) %>%
  pivot_wider(names_from = Group, values_from = proportion, values_fill = 0)

meta <- theta.collapsed %>% dplyr::select(!starts_with("Group"))
cluster <- theta.collapsed %>% dplyr::select(starts_with("Group"))

cluster <- compositions::ilr(cluster)
colnames(cluster) <- paste0("Group", 1:ncol(cluster))
data <- cbind(meta, cluster)
##### Fit Manova.RM #######
response_vars <- grep("^Group", colnames(data), value = TRUE)
fit <- multRM(as.formula(paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ Response*Time")),
              data = data, subject = "Sample", within = "Time", iter = 1000)
saveRDS(fit, file = paste0(dir, "/", thetaManova_dir, "/Manova_sims_",set_level,".rds"))
cat("Complete SIM Manova \n")
