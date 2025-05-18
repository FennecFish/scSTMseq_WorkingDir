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
library(mclust)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
design <- "nSample20_nCellType10_Batch_CancerCell/"
scSTM <- "scSTM_Pooled_Content_Sample_Prevalence_Time/"
path = paste0(dir, design, scSTM)
files <- list.files(path)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

# nCellType10_1732766297_HighVar1,nCellType10_1732766292_HighVar1, not significant
# nCellType10_1732766283_HighVar1, nCellType10_1732766268_HighVar1, significant
name <- c("scSTM_1732766297_HighVar1.rds", "scSTM_1732766292_HighVar1.rds",
          "scSTM_1732766283_HighVar1.rds", "scSTM_1732766268_HighVar1.rds")
files <- files[which(files %in% name)]
scSTMobj <- readRDS(paste0(path, files[1]))
# the following function will read in scSTMobj, and then calculate proportion difference for scSTM and simulated data
prop_diff <- function(scSTMobj, nsims = 5){
  # select top scSTMobj
  model <- select_top_scSTM(scSTMobj)
  # perform Manova.RM
  PosteriorRep <- PosteriorPropRep(model = model, nsims = nsims, Sample = "Sample", Time = "Time", Group = NULL)
  Manova.res <- ThetaManova(PosteriorRep, Time = "Time", Group = NULL)
  metadata <- colData(model$settings$sce) %>% as.data.frame() %>% dplyr::select(Cell, Sample, Time, Response)
  
  collaposed.diff <- lapply(PosteriorRep, function(x){
    # take difference between timepoints
    diff <- x %>%
      group_by(Sample) %>%
      summarise(across(starts_with("V"), ~ .[Time == "Time2"] - .[Time == "Time1"]), .groups = 'drop') %>%
      summarise(across(starts_with("V"), mean), .groups = 'drop')
    return(diff)
  })
  
  #### simulation difference
  metadata <- colData(model$settings$sce) %>% as.data.frame() 
  
  group_proportions <- metadata %>%
    group_by(Sample, Time, Response, Group) %>%
    summarise(total_sum = sum(sum), .groups = 'drop') %>%
    group_by(Sample, Time) %>%
    mutate(group_proportion = total_sum / sum(total_sum)) %>%
    ungroup()
  
  proportion_differences <- group_proportions %>%
    pivot_wider(
      names_from = Group,
      values_from = group_proportion) %>% 
    group_by(Sample, Time, Response) %>%
    summarise(
      across(starts_with("Group"), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  mean_diff <- proportion_differences %>%
    group_by(Sample, Response) %>%
    summarise(across(starts_with("Group"), ~ .[Time == "Time2"] - .[Time == "Time1"]), .groups = 'drop') %>%
    group_by(Response) %>%
    summarise(across(starts_with("Group"), mean), .groups = 'drop')
  
  # ARI
  pooled_cluster <- cluster_scSTM(model)
  metadata$Cluster <- pooled_cluster[match(rownames(metadata), names(pooled_cluster))]
  adjusted_rand_indices <- sapply(metadata %>% dplyr::select(ends_with("Cluster")), function(x) {
    adjustedRandIndex(x, metadata$Group)
  })
  return(list(ARI = adjusted_rand_indices,
              MANOVA = Manova.res,
              scSTM_diff = collaposed.diff,
              sims_diff = mean_diff))
}


x1 <- prop_diff(scSTMobj, nsims = 2)


# PosteriorRep <- PosteriorPropRep(model = model, nsims = 2, Sample = "Sample", Time = "Time", Group = "Response")
# Manova.res <- ThetaManova(PosteriorRep, Time = "Time", Group = "Response")

