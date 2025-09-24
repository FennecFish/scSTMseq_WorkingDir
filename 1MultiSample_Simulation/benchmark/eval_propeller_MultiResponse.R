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

#################################################################################
################################# ARI #########################################
#################################################################################
ari_file <- list.files("res/1MultiSample_MultiResponse_Simulation/seurat_ARI")
dat <- data.frame()
for (file in ari_file){
  temp <- readRDS(paste0("res/1MultiSample_MultiResponse_Simulation/seurat_ARI/", file))
  dat <- rbind(dat, temp)
}
colnames(dat) <- c("design", "set_level", "ARI")
res.temp <- data.frame(
  modelType = sapply(strsplit(dat$set_level, "_"), `[`, 2),
  seed = sapply(strsplit(dat$set_level, "_"), `[`, 1),
  nCellType = as.numeric(gsub("nCellType", "", sapply(strsplit(dat$design, "_"), `[`, 2))),
  nSample = as.numeric(gsub("nSample", "", sapply(strsplit(dat$design, "_"), `[`, 1))),
  Batch = ifelse(sapply(strsplit(dat$design, "_"), `[`, 3)=="noBatch", FALSE, TRUE),
  CancerType = ifelse(sapply(strsplit(dat$design, "_"), `[`, 4)=="StromalCell", FALSE, TRUE),
  ARI = as.numeric(dat$ARI)
)
write.csv(res.temp, file = "res/1MultiSample_MultiResponse_Simulation/seurat_ARI.csv")

#################################################################################
############################### Propeller #######################################
#################################################################################
dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse"
propeller.dir <- "propeller_BlockSample"
paths <- Sys.glob(file.path(dir, "*", propeller.dir)) 
files <- unlist(lapply(paths, function(path) list.files(path, full.names = TRUE)))

res <- vector(mode = "list")
for(file_name in files){
  set_level <- sub("^[^_]*_(.*)\\.rds$", "\\1",  basename(file_name))
  design <- sub(".*/MultiResponse/([^/]+)/.*", "\\1", file_name)
  nCellType = unlist(strsplit(design, "_"))[2]
  # asin <- readRDS(paste0(dir, "/", design, "/", propeller.dir, "/asin_", set_level, ".rds"))
  logit <- readRDS(paste0(dir, "/", design, "/", propeller.dir, "/logit_", set_level, ".rds"))
  temp <- logit$p.value %>%
    as.data.frame() %>%
    mutate(across(everything(), ~ p.adjust(.x, method = "BH")))
  
  list.name <- paste0(nCellType, "_", set_level)
  res[[list.name]] <- temp
  cat(file_name, "\n")
}
saveRDS(res, file = paste0("res/1MultiSample_MultiResponse_Simulation/", propeller.dir, "seurat.rds"))

#################################################################################
################################# Power #########################################
#################################################################################
dat <- readRDS("res/1MultiSample_MultiResponse_Simulation/propeller_BlockSampleseurat.rds")
threshold <- c(0.01, 0.05, 0.1)
dat.fdr <- lapply(names(dat), function(name){
  x <- dat[[name]]
  # asin.sig.change <- sapply(threshold, function(t) x$fdr.asin < t)
  # asin.sig.change <- colSums(asin.sig.change) != 0
  
  Response_Change <- sapply(threshold, function(t) x$Response_Change < t)
  Response_Change <- colSums(Response_Change) != 0
  
  nonResponse_Change <- sapply(threshold, function(t) x$nonResponse_Change < t)
  nonResponse_Change <- colSums(nonResponse_Change) != 0
  temp <- rbind(Response_Change, nonResponse_Change)
  colnames(temp) <- threshold
  return(temp)
})
names(dat.fdr) <- names(dat)

topic_plot <- do.call(rbind, lapply(names(dat.fdr), function(name) {
  df <- dat.fdr[[name]] %>% as.data.frame()
  df$Type <- rownames(df)
  df$EffectSize <- ifelse(unlist(strsplit(name, "_"))[3] == "NullModel", 0 ,as.numeric(sub(".*HighVar([0-9.]+)$", "\\1", name)))
  df$nCellType <- as.numeric(n_cell_types <- sub(".*nCellType(\\d+)_.*", "\\1", name))
  return(df)
}))

topic_plot <- topic_plot %>% 
  group_by(EffectSize, nCellType, Type) %>% 
  summarize(Power_0.01 = mean(`0.01`), 
            Power_0.05 = mean(`0.05`), 
            Power_0.1 = mean(`0.1`),
            .groups = "drop")
topic_plot <- topic_plot %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Threshold", values_to = "Power") %>%
  mutate(Threshold = as.numeric(gsub("Power_", "", Threshold)))

# png("res/1MultiSample_SingleResponse_Simulation/Power_Propeller_block_limma_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png", 
#     width = 2500, height = 1800, res = 220)
ggplot(topic_plot, aes(x = Threshold, y = Power, color = Type, group = Type)) +
  geom_point(size = 3) +
  facet_grid(nCellType~EffectSize, 
             labeller = labeller(EffectSize = function(x) paste("Effect Size =", x),
                                 nCellType = function(x) paste("nCellType =", x))) +
  labs(title = "Power for Single Response Simulation Using Propeller",
       subtitle = "Time Variable Included with Pooled Gamma Model",
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
# dev.off()

#################################################################################
############################# Type I Error ######################################
#################################################################################
# dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_limma_SingleResponse_NullHypothesis.rds")
dat <- readRDS("res/1MultiSample_SingleResponse_Simulation/1000sims/1000Propeller_block_limma_SingleResponse_NullHypothesis.rds")
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

png("res/1MultiSample_SingleResponse_Simulation/1000sims/TypeIError_1000Propeller_block_limma_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png", 
    width = 2500, height = 1800, res = 220)
ggplot(topic_plot, aes(x = Threshold, y = TypeIError, color = Type, group = Type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "red") +
  labs(title = "Type I Error for Single Response Simulation Using Propeller",
       subtitle = "Time Variable Included with Pooled Gamma Model",
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
