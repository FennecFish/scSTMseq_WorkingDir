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

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
propeller.dir <- "propeller_limma_sample_block"
paths <- Sys.glob(file.path(dir, "*nSample10*", propeller.dir)) 
files <- unlist(lapply(paths, function(path) list.files(path, pattern = "EffectSize", full.names = TRUE)))

res <- vector(mode = "list")
for(file_name in files){
  set_level <- sub("^[^_]*_(.*)\\.rds$", "\\1",  basename(file_name))
  design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
  nSample = unlist(strsplit(design, "_"))[1]
  nCellType = unlist(strsplit(design, "_"))[2]
  Batch = unlist(strsplit(design, "_"))[3]
  CancerCell = unlist(strsplit(design, "_"))[4]
  asin <- readRDS(paste0(dir, "/", design, "/", propeller.dir, "/asin_", set_level, ".rds"))
  logit <- readRDS(paste0(dir, "/", design, "/", propeller.dir, "/logit_", set_level, ".rds"))

  temp <- data.frame(Cluster = rownames(asin),
                     coefficients.asin = asin$coefficients,
                     coefficients.logit = logit$coefficients,
                     p.value.asin = asin$p.value,
                     p.value.logit = logit$p.value,
                     fdr.asin = asin$fdr,
                     fdr.logit = logit$fdr)
  list.name <- paste0(nSample, "_", nCellType, "_", Batch, "_", CancerCell, "_", set_level)
  res[[list.name]] <- temp
  cat(file_name, "\n")
}
# set <- basename(dirname(paths))[1]
# saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_",set, "_NullHypothesis.rds"))
saveRDS(res, file = paste0("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample10_EffectSize_AltHypothesis.rds"))
#################################################################################
################################# Power #########################################
#################################################################################
dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample5_NullHypothesis.rds")
dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample10_AllHypothesis.rds")
dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/Propeller_BlockSample_nSample10_EffectSize_AltHypothesis.rds")
threshold <- c(0.01, 0.05, 0.1)
# use minimum pValue
dat.fdr <- lapply(names(dat), function(name){
  x <- dat[[name]]
  
  # use fdr rate
  asin.fdr.change <- sapply(threshold, function(t) x$fdr.asin < t)
  asin.fdr.change <- colSums(asin.fdr.change) != 0

  logit.fdr.change <- sapply(threshold, function(t) x$fdr.logit < t)
  logit.fdr.change <- colSums(logit.fdr.change) != 0

  # use minimum pvalue
  asin.pmin.change <- sapply(threshold, function(t)  min(x$p.value.asin) < t)
  logit.pmin.change <- sapply(threshold, function(t)  min(x$p.value.logit) < t)
  
  temp <- rbind(asin.fdr.change, logit.fdr.change, asin.pmin.change, logit.pmin.change)
  colnames(temp) <- threshold
  return(temp)
})
names(dat.fdr) <- names(dat)

topic_plot <- do.call(rbind, lapply(names(dat.fdr), function(name) {
  df <- dat.fdr[[name]]
  fdr.asin <- df[rownames(df) == "asin.fdr.change", , drop = FALSE]
  fdr.logit <- df[rownames(df) == "logit.fdr.change", , drop = FALSE]
  pmin.asin <- df[rownames(df) == "asin.pmin.change", , drop = FALSE]
  pmin.logit <- df[rownames(df) == "logit.pmin.change", , drop = FALSE]
  temp <- as.data.frame(rbind(fdr.asin, fdr.logit, pmin.asin, pmin.logit))
  temp$Type <-rep(c("asin", "logit"), 2)
  temp$pValueType <-rep(c("FDR", "pValue.Min"), each = 2)

  if(grepl("Null", name)){
    EffectSize <- 0
  }else{
    # EffectSize <- as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", name))
    EffectSize <- as.numeric(sub(".*EffectSize(-?[0-9.]+)$", "\\1", name))
  }
  temp$EffectSize <- EffectSize
  nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", name))
  temp$nCellType <- nCellType
  
  temp$nSample = as.numeric(sub(".*nSample(\\d+)_.*", "\\1", name))
  temp$Batch = unlist(strsplit(name, "_"))[3]
  temp$CancerCell = unlist(strsplit(name, "_"))[4]

  return(temp)
}))

topic_plot <- topic_plot %>%
  group_by(EffectSize, nCellType, Type, pValueType, nSample, Batch, CancerCell) %>%
  summarize(Power_0.01 = mean(`0.01`),
            Power_0.05 = mean(`0.05`),
            Power_0.1 = mean(`0.1`),
            .groups = "drop")
topic_plot <- topic_plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "Power") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)))

topic_plot <- topic_plot %>%
  filter(Batch == "noBatch", CancerCell == "CancerCell", nCellType %in% c(8, 10, 12))
png("res/1MultiSample_SingleResponse_Simulation/Power_Propeller_BlockReplicates_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 2500, height = 1800, res = 220)
ggplot(topic_plot %>% filter(pValueType == "FDR"), aes(x = Threshold, y = Power, color = Type, group = Type)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize,
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Power for Single Response Simulation Using Propeller",
       subtitle = "limma with block replicates",
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
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_limma_SingleResponse_NullHypothesis.rds")
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_block_limma_SingleResponse_NullHypothesis.rds")
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_BlockReplicates_SingleResponse_NullHypothesis.rds")
threshold <- c(0.01, 0.05, 0.1)
dat.fdr <- lapply(names(dat), function(name){
  x <- dat[[name]]
  asin.sig.change <- sapply(threshold, function(t) x$fdr.asin < t)
  asin.sig.change <- colSums(asin.sig.change) != 0

  logit.sig.change <- sapply(threshold, function(t) x$fdr.logit < t)
  logit.sig.change <- colSums(logit.sig.change) != 0

  temp <- rbind(asin.sig.change, logit.sig.change)
  colnames(temp) <- threshold
  return(temp)
})
names(dat.fdr) <- names(dat)

topic_plot <- do.call(rbind, lapply(names(dat.fdr), function(name) {
  df <- dat.fdr[[name]]
  df.asin <- df[rownames(df) == "asin.sig.change", , drop = FALSE]
  df.logit <- df[rownames(df) == "logit.sig.change", , drop = FALSE]
  temp <- as.data.frame(rbind(df.asin, df.logit))
  temp$Type <- c("asin", "logit")
  temp$EffectSize <- "Null"
  nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", name))
  temp$nCellType <- nCellType
  return(temp)
}))

topic_plot <- topic_plot %>%
  group_by(nCellType, Type) %>%
  summarize(Power_0.01 = mean(`0.01`),
            Power_0.05 = mean(`0.05`),
            Power_0.1 = mean(`0.1`),
            .groups = "drop")
topic_plot <- topic_plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "TypeIError") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)))

png("res/1MultiSample_SingleResponse_Simulation/1000sims/TypeIError_1000Propeller_BlockReplicates_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
    width = 2500, height = 1800, res = 220)
ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Type, group = Type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Type I Error for Single Response Simulation Using Propeller",
       subtitle = "limma with block replicates",
       x = "Thresholds",
       y = "Type I Error") +
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
