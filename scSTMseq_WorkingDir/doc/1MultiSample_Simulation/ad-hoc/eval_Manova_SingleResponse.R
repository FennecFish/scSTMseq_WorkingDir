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
# manova.dir <- "1000ManovaTheta_Pooled_noContent_Prevalence_Time"
manova.dir <- "ManovaRM"
paths <- Sys.glob(file.path(dir, "*nSample10*noBatch_CancerCell*", manova.dir)) 
# files <- unlist(lapply(paths, function(path) list.files(path, pattern = "Manova_pValue_.*EffectSize-?[0-9.]+\\.rds", full.names = TRUE)))
files <- unlist(lapply(paths, function(path) list.files(path, pattern = "Manova_pValue_", full.names = TRUE)))

res <- vector(mode = "list")

for(file_name in files){
  set_level <- sub("^Manova_pValue_(.*)\\.rds$", "\\1", basename(file_name))
  design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
  nSample = unlist(strsplit(design, "_"))[1]
  nCellType = unlist(strsplit(design, "_"))[2]
  Batch = unlist(strsplit(design, "_"))[3]
  CancerCell = unlist(strsplit(design, "_"))[4]
  
  
  # Attempt to read the .rds file, skipping on error
  temp <- tryCatch({
    ThetaManova <- readRDS(file_name)
    data.frame(pValue.mats = ThetaManova$rm.mats)
  }, error = function(e) {
    cat("Error reading file:", file_name, "- Skipping\n")
    return(NULL)
  })
  
  # If readRDS failed, skip to the next iteration
  if (is.null(temp)) {
    next
  }
  
 #  temp <- data.frame(pValue.mats = ThetaManova$rm.mats)
  colnames(temp) <- c("pValue.mats")
  # colnames(temp) <- c("pValue.wts", "pValue.mats")
  list.name <- paste0(nSample, "_", nCellType, "_", Batch, "_", CancerCell, "_", set_level)
  res[[list.name]] <- temp
  cat(file_name, "\n")
}
# name <- sub("^ManovaTheta_", "",basename(manova.dir))
# set <- basename(dirname(paths))[1]
saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample10_AllHypothesis.rds"))

# #################################################################################
# ################################## plot #########################################
# #################################################################################
# plot_data_generate <- function(thresholds, pValueList, by_group = NULL){
#   exist.change.list <- lapply(thresholds, function(thr) {
#     res <- pValueList %>%
#       rownames_to_column("Model") %>%
#       mutate(nSample = as.numeric(sub(".*nSample(\\d+)_.*", "\\1", Model)),
#              nCellType = as.numeric(sub(".*nCellType(\\d+)_.*", "\\1", Model)),
#              # EffectSize = ifelse(grepl("Null", Model), 0, as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", Model))),
#              EffectSize = as.numeric(sub(".*EffectSize(-?[0-9.]+)$", "\\1", Model)),
#              Batch = unlist(strsplit(Model, "_"))[3],
#              CancerCell = unlist(strsplit(Model, "_"))[4],
#              wts.sig = wts < thr,
#              mats.sig = mats < thr)
#     return(res)
#   })
# 
#   ProportionOfTrue <- data.frame()
#   for (i in seq_along(thresholds)) {
#     threshold <- thresholds[i]
#     if(is.null(by_group)){
# 
#       prop_truth <- exist.change.list[[i]] %>%
#         as.data.frame() %>%
#         group_by(EffectSize, nCellType, nSample, Batch, CancerCell) %>%
#         summarize(
#           wts_sig_proportion = mean(wts.sig),
#           mats_sig_proportion = mean(mats.sig),
#           .groups = "drop"
#         ) %>%
#         mutate(Threshold = threshold)
# 
# 
#     }else{
#       # if we want to group the power/ type I error by effect size
#       prop_truth <- exist.change.list[[i]] %>%
#         as.data.frame() %>%
#         group_by(EffectSize, nCellType, nSample, Batch, CancerCell) %>%
#         summarize(
#           wts_sig_proportion = mean(wts.sig),
#           mats_sig_proportion = mean(mats.sig),
#           .groups = "drop"
#         ) %>%
#         mutate(Threshold = threshold)
#     }
#     # Store the results in a dataframe
#     ProportionOfTrue <- bind_rows(ProportionOfTrue, prop_truth)
#   }
#   return(ProportionOfTrue)
# }
# 
# #################################################################################
# ######################### Type I Error & Power ##################################
# #################################################################################
# res <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample5.rds")
# res <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample10.rds")
# res <- readRDS("res/1MultiSample_SingleResponse_Simulation/Manova_pValue_nSample10_EffectSize.rds")
# manova.pValue <- lapply(res, function(x) {
#   data.frame(wts =quantile(x$pValue.wts, probs = 0.50),
#              mats = quantile(x$pValue.mats, probs = 0.50))
# })
# manova.pValue <- do.call(rbind, manova.pValue)
# 
# threshold <- c(0.01, 0.05, 0.1)
# topic_plot <- plot_data_generate(threshold, pValueList = manova.pValue)
# topic_plot <- topic_plot %>%
#   dplyr::rename(WTS = wts_sig_proportion,
#          MATS = mats_sig_proportion) %>%
#   pivot_longer(cols = c("WTS", "MATS"), names_to = "Methods", values_to = "TypeIError")
# 
# topic_plot <- topic_plot %>%
#   filter(Batch == "noBatch", CancerCell == "CancerCell", Methods == "MATS")
# png("res/1MultiSample_SingleResponse_Simulation/TypeIError_SingleResponse_Batch_CancerCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 1500, height = 1200, res = 220)
# ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Methods, group = Methods)) +
#   facet_grid(nCellType~EffectSize) +
#   geom_point(size = 3) +
#   geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
#   labs(title = "Type I Error for Single Response Simulation Under Various Thresholds",
#        subtitle = "Time Covariate Included, No Batch Effect and Stromal Cells Only",
#        x = "Threshold",
#        y = "Type I Error") +
#   theme_bw() +
#   scale_x_continuous(breaks = threshold) +  # Ensure proper thresholds on X axis
#   theme(
#     plot.title = element_text(size = 15, hjust = 0.5),  # Increase title text size
#     axis.title.x = element_text(size = 14),  # Increase X axis title size
#     axis.title.y = element_text(size = 14),  # Increase Y axis title size
#     axis.text.x = element_text(size = 14),   # Increase X axis labels size
#     axis.text.y = element_text(size = 14)   # Increase Y axis labels size
#   )
# dev.off()
# 
# png("res/1MultiSample_SingleResponse_Simulation/Power_SingleResponse_Batch_CancerCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 2500, height = 1800, res = 220)
# ggplot(topic_plot, aes(x = Threshold, y = Power, color = Methods, group = Methods)) +
#   geom_point(size = 3) +
#   facet_grid(~EffectSize,
#              labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
#                                  nCellType = function(x) paste("nCellType =", x))) +
#   labs(title = "Power for Single Response Simulation with Batch Effect and Cancer-like Cells",
#        subtitle = "Summing Across Theta",
#        x = "Thresholds",
#        y = "Power") +
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
