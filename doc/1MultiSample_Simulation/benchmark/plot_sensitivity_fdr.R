setwd("/proj/milovelab/wu/scSTMseq_WorkingDir")
library(Matrix)
library(dplyr)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)
library(MASS)
source("doc/1MultiSample_Simulation/benchmark/utils_functions.R")

threshold <- c(0.01, 0.05, 0.1)

scSTMseq <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample10_AllHypothesis.rds")
scSTMseq <- scSTMseq[grep("noBatch_CancerCell",names(scSTMseq))]
# take the medium
scSTMseq.pValue <- do.call(rbind,lapply(scSTMseq, function(x) {
  data.frame(mats = quantile(x$pValue.mats, probs = 0.50))
}))
scSTMseq_plot <- scSTMseq_plot_data(threshold, pValueList = scSTMseq.pValue)

propeller.sample.block <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample10_AllHypothesis.rds")
propeller.sample.block <- propeller.sample.block[grep("noBatch_CancerCell",names(propeller.sample.block))]
propeller.sample.block.plot <- propeller_plot_dat(dat = propeller.sample.block, threshold = threshold)

plot_data <- bind_rows(scSTMseq_plot, propeller.sample.block.plot) %>%
  filter(nCellType %in% c(8, 10, 12))

png("res/1MultiSample_SingleResponse_Simulation/FDRvsSensivity_PropellerVSscSTMseq_nSample10_noBatch_CancerCell.png",
    width = 3000, height = 1800, res = 220)
ggplot(plot_data, aes(x = FDR, y = Power, color = Method)) +
  facet_wrap( ~ nCellType * EffectSize,
    nrow = 3,
    labeller = labeller(
      nCellType = function(x) paste("nCellType =", x),
      EffectSize = function(x) paste("Effect Size of Change =", x))) +
  geom_point(cex = 2) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  # geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(
    x = "FDR",
    y = "Sensitivity",
    title = "FDR vs. Sensitivity Across Different Cell Types and Effect Size"
  ) +
  theme_bw()
dev.off()


scSTMseq <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_SingleResponse_noBatch_StromalCell_Pooled_noContent_Prevalence_Time.rds")
scSTMseq <- scSTMseq[grep("HighVar",names(scSTMseq))]
propeller <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_limma_SingleResponse_AltHypothesis.rds")
propeller <- propeller[grep("HighVar",names(propeller))]
propeller.replicates.block <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockReplicates_SingleResponse_AltHypothesis.rds")
propeller.replicates.block <- propeller.replicates.block[grep("HighVar",names(propeller.replicates.block))]
propeller.sample.block <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_block_limma_SingleResponse_AltHypothesis.rds")
propeller.sample.block <- propeller.sample.block[grep("HighVar",names(propeller.sample.block))]

threshold <- c(0.01, 0.05, 0.1)
### manipulate scSTMseq
scSTMseq.pValue <- lapply(scSTMseq, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
scSTMseq.pValue <- do.call(rbind, scSTMseq.pValue)
scSTMseq.pValue$EffectSize <- as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", rownames(scSTMseq.pValue)))
scSTMseq.pValue$nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", rownames(scSTMseq.pValue)))
scSTMseq_plot <- plot_data_generate(threshold, pValueList = scSTMseq.pValue, by_group = "EffectSize")
scSTMseq_plot <- scSTMseq_plot %>%
  dplyr::rename(WTS = wts_sig_proportion,
                MATS = mats_sig_proportion) %>%
  pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "Power") %>%
  mutate(Methods = paste0("scSTMseq.", Methods)) %>%
  dplyr::select(EffectSize, nCellType, Threshold, Power, Methods)

# data clean for propeller
propeller.plot <- propeller_plot_dat(dat = propeller, threshold = threshold)
propeller.plot <- propeller.plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "Power") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
         Methods = paste0("propeller.unpaired.",Type)) %>%
  dplyr::select(-Type) 

propeller.replicates.block.plot <- propeller_plot_dat(dat = propeller.replicates.block, threshold = threshold)
propeller.replicates.block.plot <- propeller.replicates.block.plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "Power") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
         Methods = paste0("propeller.BlockRep.",Type)) %>%
  dplyr::select(-Type)

propeller.sample.block.plot <- propeller_plot_dat(dat = propeller.sample.block, threshold = threshold)
propeller.sample.block.plot <- propeller.sample.block.plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "Power") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
         Methods = paste0("propeller.BlockSample.",Type)) %>%
  dplyr::select(-Type)

dat <- rbind(scSTMseq_plot, propeller.plot, propeller.replicates.block.plot, propeller.sample.block.plot)

method.sub <- unique(dat$Methods)[c(2,4,6,8)]
png("res/1MultiSample_SingleResponse_Simulation/Power_PropellerVSscSTMseq_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 3000, height = 1800, res = 220)
ggplot(dat %>% filter(Methods %in% method.sub), aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize,
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Group Simulation Comparing scSTMseq and Propeller",
       # subtitle = "limma with block replicates",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14),    # Increase facet label text size
    legend.position = "bottom",             # Move legend to the bottom
    legend.text = element_text(size = 14),  # Increase legend text size
    panel.spacing.x = unit(1, "lines"),   # Horizontal spacing (unchanged, default small spacing)
    panel.spacing.y = unit(0.5, "lines")    # Vertical spacing (increased)
  )
dev.off()

png("res/1MultiSample_SingleResponse_Simulation/Power_AllMethods_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 3000, height = 1800, res = 220)
ggplot(dat, aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize,
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Group Simulation Across Methods",
       # subtitle = "limma with block replicates",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14),    # Increase facet label text size
    legend.position = "bottom",             # Move legend to the bottom
    legend.text = element_text(size = 14),  # Increase legend text size
    panel.spacing.x = unit(1, "lines"),   # Horizontal spacing (unchanged, default small spacing)
    panel.spacing.y = unit(0.5, "lines")    # Vertical spacing (increased)
  )
dev.off()

method.sub <- unique(dat$Methods)[3:6]
png("res/1MultiSample_SingleResponse_Simulation/Power_Propeller_RepVSnoRep_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 2500, height = 1800, res = 220)
ggplot(dat %>% filter(Methods %in% method.sub ), aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize,
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Response Simulation Using Propeller",
       subtitle = "Comparing Inclusion of Clustering Uncertainty",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
dev.off()

method.sub <- unique(dat$Methods)[c(2, 6, 8)]
png("res/1MultiSample_SingleResponse_Simulation/Power_PropellerVSscSTMseq_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 2500, height = 1800, res = 220)
ggplot(dat %>% filter(Methods %in% method.sub), aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize,
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Response Simulation Using Propeller vs scSTMseq",
       # subtitle = "Comparing Inclusion of Clustering Uncertainty",
       x = "Thresholds",
       y = "Power") +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
    axis.title.x = element_text(size = 16),  # Increase X axis title size
    axis.title.y = element_text(size = 16),  # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),   # Increase Y axis labels size
    strip.text = element_text(size = 14)  # Increase facet label text size
  )
dev.off()
#################################################################################
############################# Type I Error ######################################
#################################################################################
scSTMseq <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample10_AllHypothesis.rds")
scSTMseq <- scSTMseq[grep("NullModel",names(scSTMseq))]
scSTMseq <- scSTMseq[grep("noBatch_CancerCell",names(scSTMseq))]
# propeller <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_limma_SingleResponse_NullHypothesis.rds")
# propeller <- propeller[grep("NullModel",names(propeller))]
propeller.sample.block <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample10_AllHypothesis.rds")
propeller.sample.block <- propeller.sample.block[grep("NullModel",names(propeller.sample.block))]
propeller.sample.block <- propeller.sample.block[grep("noBatch_CancerCell",names(propeller.sample.block))]
# propeller.replicates.block <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_BlockReplicates_SingleResponse_NullHypothesis.rds")
# propeller.replicates.block <- propeller.replicates.block[grep("NullModel",names(propeller.replicates.block))]
threshold <- c(0.01, 0.05, 0.1)

### manipulate scSTMseq
scSTMseq.pValue <- lapply(scSTMseq, function(x) {
  data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
             mats = quantile(x$pValue.mats, probs = 0.50))
})
scSTMseq.pValue <- do.call(rbind, scSTMseq.pValue)
scSTMseq_plot <- plot_data_generate(threshold, pValueList = scSTMseq.pValue, by_group = NULL)
scSTMseq_plot <- scSTMseq_plot %>%
  dplyr::rename(WTS = wts_sig_proportion,
                MATS = mats_sig_proportion) %>%
  pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "TypeIError") %>%
  mutate(Methods = paste0("scSTMseq.", Methods)) %>%
  dplyr::select(Threshold, nCellType, TypeIError, Methods)

# # data clean for propeller
# propeller.plot <- propeller_plot_dat(dat = propeller, threshold = threshold)
# propeller.plot <- propeller.plot %>%
#   pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "TypeIError") %>%
#   mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
#          Methods = paste0("propeller.unpaired.",Type)) %>%
#   dplyr::select(Threshold, TypeIError, Methods)
# 
# propeller.replicates.block.plot <- propeller_plot_dat(dat = propeller.replicates.block, threshold = threshold)
# propeller.replicates.block.plot <- propeller.replicates.block.plot %>%
#   pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "TypeIError") %>%
#   mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
#          Methods = paste0("propeller.BlockRep.",Type)) %>%
#   dplyr::select(Threshold, TypeIError, Methods)

propeller.sample.block.plot <- propeller_plot_dat(dat = propeller.sample.block, threshold = threshold)
propeller.sample.block.plot <- propeller.sample.block.plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "TypeIError") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)),
         Methods = paste0("propeller.BlockSample.",Type)) %>%
  dplyr::select(Threshold, nCellType, TypeIError, Methods)

# dat <- rbind(scSTMseq_plot, propeller.plot, propeller.replicates.block.plot, propeller.sample.block.plot)
dat <- rbind(scSTMseq_plot, propeller.sample.block.plot)

method.sub <- unique(dat$Methods)[c(2,3,4)]

png("res/1MultiSample_SingleResponse_Simulation//TypeIError_PropellerVSscSTMseq_nSample10_noBatch_CancerCell.png",
    width = 2500, height = 1800, res = 220)
ggplot(dat %>% filter(Methods %in% method.sub, nCellType %in% c(8, 10, 12)), 
       aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
  facet_grid(~nCellType, labeller = labeller(nCellType = function(x) paste0("nCellType = ", x))) +
  geom_point(size = 3) +
  # geom_jitter(size = 3, width = 0.002, height = 0, alpha = 1) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(
    title = "Type I Error Using scSTMseq and Propeller with Paired Sample",
    x = "Threshold",
    y = "Type I Error",
    color = "Method"  # Renaming the legend title
  ) +
  theme_bw() +
  scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on the X axis
  scale_color_manual(
    values = c("scSTMseq.MATS" = "skyblue", "propeller.BlockSample.logit" = "orange",
               "propeller.BlockSample.asin" = "darkgrey"),  # Set custom colors
    labels = c("scSTMseq.MATS" = "scSTMseq", "propeller.BlockSample.logit" = "Propeller.logit",
               "propeller.BlockSample.asin" = "Propeller.asin")  # Rename legend entries
  ) +
  theme(
    plot.title = element_text(size = 18, hjust = 0),  # Left-align title
    axis.title.x = element_text(size = 16),          # Increase X axis title size
    axis.title.y = element_text(size = 16),          # Increase Y axis title size
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),           # Increase Y axis labels size
    strip.text = element_text(size = 14),            # Increase facet label text size
    legend.position = "bottom",                      # Move legend to the bottom
    legend.text = element_text(size = 14)
  )
dev.off()
# 
# png("res/1MultiSample_SingleResponse_Simulation/1000sims/TypeIError_AllMethods_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 2500, height = 1800, res = 220)
# ggplot(dat, aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
#   geom_point(size = 3) +
#   geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
#   labs(title = "Type I Error for Single Group Simulation for All Methods",
#        #subtitle = "Time Covariate Included, No Batch Effect and Stromal Cells Only",
#        x = "Threshold",
#        y = "Type I Error") +
#   theme_bw() +
#   scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
#   theme(
#     plot.title = element_text(size = 18, hjust = 0),  # Left-align title
#     plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
#     axis.title.x = element_text(size = 16),  # Increase X axis title size
#     axis.title.y = element_text(size = 16),  # Increase Y axis title size
#     axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
#     axis.text.y = element_text(size = 14),   # Increase Y axis labels size
#     strip.text = element_text(size = 14),  # Increase facet label text size
#     legend.position = "bottom",             # Move legend to the bottom
#     legend.text = element_text(size = 14)
#   )
# dev.off()
# 
# 
# 
# method.sub <- unique(dat$Methods)[c(2, 6, 8)]
# png("res/1MultiSample_SingleResponse_Simulation/1000sims/TypeIError_PropellerVSscSTMseq_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 2500, height = 1800, res = 220)
# ggplot(dat %>% filter(Methods %in% method.sub), aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
#   geom_point(size = 3) +
#   # geom_jitter(size = 3, height = 0, alpha = 1) +
#   geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
#   labs(title = "Type I Error for Single Response Simulation Using Propeller vs scSTMseq",
#        # subtitle = "Comparing Inclusion of Clustering Uncertainty",
#        x = "Threshold",
#        y = "Type I Error") +
#   theme_bw() +
#   scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
#   theme(
#     plot.title = element_text(size = 18, hjust = 0),  # Left-align title
#     plot.subtitle = element_text(size = 14, hjust = 0),  # Left-align subtitle
#     axis.title.x = element_text(size = 16),  # Increase X axis title size
#     axis.title.y = element_text(size = 16),  # Increase Y axis title size
#     axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
#     axis.text.y = element_text(size = 14),   # Increase Y axis labels size
#     strip.text = element_text(size = 14)  # Increase facet label text size
#   )
# dev.off()
