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

##### analysis ####'
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
manova.dir <- "ManovaTheta_Pooled_Content_Sample_Prevalence_Time"
paths <- Sys.glob(file.path(dir, "*nCellType10_Batch_*", manova.dir)) 
files <- unlist(lapply(paths, function(path) list.files(path, pattern = "Manova_sims", full.names = TRUE)))

res <- vector(mode = "list")

for(file_name in files){
  set_level <- sub("^Manova_sims_(.*)\\.rds$", "\\1", basename(file_name))
  design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
  nCellType = unlist(strsplit(design, "_"))[2]
  ThetaManova <- readRDS(file_name)
  ThetaManova <- ThetaManova$resampling[,2]
  list.name <- paste0(nCellType, "_", set_level)
  res[[list.name]] <- ThetaManova
  cat(file_name, "\n")
}
name <- sub("^ManovaTheta_", "",basename(manova.dir))
saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/Manova_sims_nSample20_", name,".rds"))

#################################################################################
################################## plot #########################################
#################################################################################
plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
  exist.change.list <- lapply(thresholds, function(thr) {
    res <- pValueList %>%
      mutate(time.sig = Time < thr)
             
    return(res)
  })
  ProportionOfTrue <- data.frame()

  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    if(is.null(by_group)){
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        summarize(
          wts_sig_proportion = mean(wts.sig),
          mats_sig_proportion = mean(mats.sig)
        ) %>%
        mutate(Threshold = threshold)
    }else{
      # if we want to group the power/ type I error by effect size
      prop_truth <- exist.change.list[[i]] %>%
        as.data.frame() %>%
        group_by(EffectSize, nCellType) %>%
        summarize(
          time_sig_proportion = mean(time.sig),
          .groups = "drop"
        ) %>%
        mutate(Threshold = threshold)
    }
    # Store the results in a dataframe
    ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
  }
  return(ProportionOfTrue)
}

dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_sims_nSample20_Pooled_Content_Sample_Prevalence_Time.rds")
dat <- do.call(rbind, dat)
colnames(dat) <- "Time"
dat <- dat %>%
  as.data.frame() %>%
  mutate(
    nCellType = as.numeric(sub("nCellType", "", sapply(strsplit(rownames(dat), "_"), `[`, 1))),
    EffectSize = sapply(strsplit(rownames(dat), "_"), `[`, 3))
dat <- dat %>% 
  mutate(EffectSize = case_when(
      EffectSize == "NullModel" ~ 0,
      grepl("HighVar", EffectSize) ~ as.numeric(sub("HighVar", "", EffectSize))
    ))

threshold <- c(0.01, 0.05, 0.1)
topic_plot <- plot_data_generate(threshold, pValueList = dat, by_group = "EffectSize")
topic_plot <- topic_plot %>%
  mutate(Threshold = as.factor(Threshold)) %>%
  dplyr::rename(Time = time_sig_proportion) %>%
  pivot_longer(cols = c("Time"), names_to = "Covariate", values_to = "Power")

# png("res/1MultiSample_SingleResponse_Simulation/Power_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png", 
#     width = 2500, height = 1800, res = 220)
ggplot(topic_plot, aes(x = Covariate, y = Power, color = Threshold, group = Threshold)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize, 
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Response Simulation with Batch Effect and Cancer-like cells",
       subtitle = "Use Simulated Data Prior to Clustering",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  # scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
# dev.off()
