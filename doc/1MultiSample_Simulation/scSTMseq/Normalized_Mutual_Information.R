setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(aricode)

# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
# design <- "nSample20_nCellType5_noBatch_StromalCell/"
scSTMdir <- "scSTM_Pooled_noContent_Prevalence_Time"
# lm <-"scSTM_LinearMixed_noContent_Prevalence_Time/"
# extract scSTM files from different design
paths <- Sys.glob(file.path(dir, "*", scSTMdir)) 
files <- unlist(lapply(paths, list.files, full.names = TRUE))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

res <- data.frame()

for(file_name in files){
  
  # extract the seed and effect size
  set_level <- sub("^scSTM_(.*)\\.rds$", "\\1", basename(file_name))
  # extract nSample and nCEll TYpe
  design <- sub(".*/SingleResponse/([^/]+)/.*", "\\1", file_name)
  # first read in pooled data
  pooled_scSTMobj <- safe_readRDS(file_name)
  pooled_scSTMobj <- select_top_scSTM(pooled_scSTMobj)
  sce <- pooled_scSTMobj$settings$sce
  x <- pooled_scSTMobj$settings$sce$Group
  theta <- pooled_scSTMobj$theta
  rownames(theta) <- pooled_scSTMobj$DocName
  colnames(theta) <- paste0("Topic", 1:ncol(theta))
  K <- pooled_scSTMobj$settings$dim$K
  # calcluate entropy
  collapsed.theta <- colSums(theta)
  norm.collaposed.theta <- collapsed.theta/sum(collapsed.theta)
  H.Y <- - sum(norm.collaposed.theta * log2(norm.collaposed.theta))
  # calculate H(Y|C)
  H.Y.C <- vector(mode = "list")
  for(group in levels(x)){
    # extract cells that in the given group
    cells <- sce$Cell[sce$Group == group]
    subset.theta <- theta[rownames(theta) %in% cells,]
    collapsed.sub.theta <- colSums(subset.theta)
    norm.sub.theta <- collapsed.sub.theta/sum(collapsed.sub.theta)
    H.Y.C[[group]] <- - sum(norm.sub.theta * log2(norm.sub.theta))/K
  }
  H.Y.C <- do.call(sum, H.Y.C)
  I.Y.C <- H.Y - H.Y.C
  
  C.prop <- as.vector(prop.table(table(x)))
  H.C <- - sum(C.prop*log2(C.prop))
  NMI <- 2*I.Y.C/(H.Y + H.C)
  
  dat <- colData(pooled_scSTMobj$settings$sce) %>% data.frame() 
  dat$pooled_scSTMseq <- pooled_cluster[match(rownames(dat), names(pooled_cluster))]
  
  # then read in linear mixed model
  # lm_scSTMobj_path <- paste0(dir, design, lm, "scSTM_", set_level, ".rds")
  # if(file.exists(lm_scSTMobj_path)){
  #   lm_scSTMobj <- readRDS(lm_scSTMobj_path)
  #   lm_cluster <- process_scSTM(lm_scSTMobj)
  #   dat$lm_scSTMseq <- lm_cluster[match(rownames(dat), names(lm_cluster))]
  # } else {dat$lm_scSTMseq <- NA}
  
  adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("scSTMseq")), function(x) {
    adjustedRandIndex(x, dat$Group)
  })
  
  # Create a data frame to store results
  res.temp <- data.frame(
    scSTMTyoe = scSTMdir,
    modelType = unlist(strsplit(set_level, "_"))[2],
    seed = unlist(strsplit(set_level, "_"))[1],
    nCellType = as.numeric(gsub("nCellType", "", unlist(strsplit(design, "_"))[2])),
    nSample = as.numeric(gsub("nSample", "", unlist(strsplit(design, "_"))[1])),
    Batch = ifelse(unlist(strsplit(design, "_"))[3]=="noBatch", FALSE, TRUE),
    CancerType = ifelse(unlist(strsplit(design, "_"))[4]=="StromalCell/", FALSE, TRUE),
    t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  )
  
  res <- bind_rows(res, res.temp)
  cat(file_name, "\n")
  rm(pooled_scSTMobj)
  rm(dat)
}

# pre_res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final.csv")
# res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final_update.csv")
# 
# res.update <- rbind(pre_res, res)
# res.update <- res %>%
#  # dplyr::select(sim, -matches("cluster$"), matches("cluster$")) %>%
#   dplyr::full_join(pre_res, by = c("sim", "seed", "control", "level", "nCellType", "nSample")) %>%
#  dplyr::select(sim, seed, control, level, nCellType, nSample,
#                -matches("cluster$"), matches("cluster$"))

save_path <- paste0("res/1MultiSample_SingleResponse_Simulation/adjRandIndex_",basename(dir), "_", basename(scSTMdir),".csv")
write.csv(res, file = save_path)

##################################################################################
########################### plot #################################################
##################################################################################
dat <- read.csv("res/1MultiSample_SingleResponse_Simulation/adjRandIndex_SingleResponse_scSTM_Pooled_noContent_Prevalence_Time.csv") 

dat_long <- dat %>% 
  mutate(nCellType = recode(nCellType, "5" = "5 Cell Types", "10" = "10 Cell Types", 
                            "15" = "15 Cell Types")) %>%
  gather(method, adjR, pooled_scSTMseq, factor_key=TRUE) %>%
  drop_na()  %>%
  group_by(nCellType, nSample, modelType, method) %>%
  summarise(mean_adjR = mean(adjR, na.rm = TRUE), .groups = "drop") 
# mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))

dat_long$nCellType <- factor(dat_long$nCellType, levels = c("5 Cell Types", "10 Cell Types", "15 Cell Types"))
dat_long$modelType <- as.numeric(sub("HighVar", "", dat_long$modelType))
# dat_long$nSample <- factor(dat_long$nSample, levels = c("nsample3", "nsample6", "nsample12"))

# png("res/sim_cluster_benchmark/adjustedRandIndex_lineplot_multiple_sample_benchmark_V4_figure.png", width = 1200, height = 1000, res = 150)
ggplot(dat_long, aes(x = modelType, y = mean_adjR, color = nCellType, group = nCellType)) +
  # geom_line(aes(linetype = line_type), linewidth = 1.2) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Effect Size",
    y = "Mean Adjusted Rand Index",
    title = "Benchmarking Clustering Accuracy Using Mean Adjusted Rand Index",
    subtitle = "Across Different Number of Cell Types"
  ) +
  theme_minimal(base_size = 12) +
  theme_bw() + 
  ylim(0, 1) + 
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) 
# dev.off()

png("res/sim_cluster_benchmark/adjustedRandIndex_boxplot_multiple_sample_benchmark_V4.png", width = 1200, height = 800, res = 160)
ggplot(dat_long, aes(x = method, y = mean_adjR)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # Remove outliers from the boxplot, adjust box width
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = method)) +  # Add jittered points for data distribution
  facet_grid(nSample ~ nCellType) +  # Use facet_grid for a more structured layout
  theme_minimal(base_size = 14) +  # Set a base font size for better readability
  labs(
    title = "Boxplot of Clustering Accuracy by Method Simulated Using Splatter",
    x = "Method",
    y = "Mean Adjusted Rand Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center and bold the title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust x-axis text for better readability
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title = element_text(size = 10),  # Increase axis title size
    strip.text = element_text(size = 10, face = "bold"),  # Bold facet labels
    legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
  )
dev.off()