# Goal: This script is to calculate ARI for scSTMseq 
# Author: Euphy Wu
# Note: Recommend submitting with bash script, with no input

setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(cluster)
library(tidyr)

dir <- "/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse"
# design <- "nSample20_nCellType5_noBatch_StromalCell/"
# scSTMdir <- "scSTM_Pooled_Content_Sample_Prevalence_Time"
# scSTMdir <- "scSTM_Pooled_noContent_Prevalence_Time"
scSTMdir <- "scSTM_Sensitivity_Pooled_Content_Prevalence_Time"
# extract scSTM files from different design
paths <- Sys.glob(file.path(dir, "nSample10_nCellType10_Batch_CancerCell", scSTMdir)) 
files <- unlist(lapply(paths, list.files, full.names = TRUE, pattern = "ARI"))

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
  
  # in the following code, we plot the convergence plot for each of those run
  
  run <- pooled_scSTMobj$runout
  obj <- vector(mode = "list")
  ari <- vector(mode = "list")
  method <- vector(mode = "list")
  for (i in 1:length(run)){
    obj[[i]] <- run[[i]]$convergence$bound
    method[[i]] <- run[[i]]$settings$init$mode
    ari[[i]] <- unlist(run[[i]]$ARI)
  }
  
  names(obj) <- unlist(method)
  
  dat <- bind_rows(
    lapply(seq_along(obj), function(index) {
      data.frame(
        index = index,
        method = names(obj)[[index]],
        iter = seq_along(obj[[index]]),
        value = obj[[index]],
        ari = ari[[index]]
      )
    })
  )

  ggplot(dat %>% filter(iter >= 2), aes(x = iter, y = value, color = factor(index))) +
    geom_line() +
    facet_grid(method~., scales = "free_y") +
    theme_minimal() +
    labs(x = "iterations", y = "Bound", title = "Convergence Plot by Initialization Method") +
    theme(legend.position = "none")
  # pooled_scSTMobj <- tryCatch({
  #   select_top_scSTM(pooled_scSTMobj)
  # }, error = function(e) {
  #   return(NULL)  # Don't use `next` here
  # })
  # if (is.null(pooled_scSTMobj)) {
  #   next  # Now `next` is at the correct scope — inside the loop
  # }
  for (i in 1:length(pooled_scSTMobj$runout)) {
    scSTMobj <- pooled_scSTMobj$runout[[i]]
    pooled_cluster <- cluster_scSTM(scSTMobj)
    
    dat <- colData(scSTMobj$settings$sce) %>% data.frame() 
    dat$pooled_scSTMseq <- pooled_cluster[match(rownames(dat), names(pooled_cluster))]
    adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("scSTMseq")), function(x) {
      adjustedRandIndex(x, dat$Group)
    })
    
    res.temp <- data.frame(
      scSTMType = scSTMdir,
      # modelType = unlist(strsplit(set_level, "_"))[2],
      # seed = unlist(strsplit(set_level, "_"))[1],
      modelType = unlist(strsplit(set_level, "_"))[4],
      seed = unlist(strsplit(set_level, "_"))[3],
      nCellType = as.numeric(gsub("nCellType", "", unlist(strsplit(design, "_"))[2])),
      nSample = as.numeric(gsub("nSample", "", unlist(strsplit(design, "_"))[1])),
      Batch = ifelse(unlist(strsplit(design, "_"))[3]=="noBatch", FALSE, TRUE),
      CancerType = ifelse(unlist(strsplit(design, "_"))[4]=="StromalCell", FALSE, TRUE),
      Num_iter = as.numeric(gsub("iter","",unlist(strsplit(set_level, "_"))[1])),
      Num_init = as.numeric(gsub("init","",unlist(strsplit(set_level, "_"))[2])),
      method = scSTMobj$settings$init$mode,
      init_llh = scSTMobj$convergence$bound[1],
      llh5 = scSTMobj$convergence$bound[5],
      llh10 = scSTMobj$convergence$bound[10],
      llh15 = scSTMobj$convergence$bound[15],
      llh30 = scSTMobj$convergence$bound[30],
      llh50 = scSTMobj$convergence$bound[50],
      max_llh = tail(scSTMobj$convergence$bound, n=1),
      scSTMseq = t(adjusted_rand_indices)  # Transpose to match the original data frame structure
      # Seurat = seurat_ari
    )
    
    res <- bind_rows(res, res.temp)
    # cat(file_name, "\n")
    rm(scSTMobj)
    rm(dat)
  }
  cat(file_name, "\n")
  # pooled_cluster <- cluster_scSTM(pooled_scSTMobj)
  # 
  # dat <- colData(pooled_scSTMobj$settings$sce) %>% data.frame() 
  # dat$pooled_scSTMseq <- pooled_cluster[match(rownames(dat), names(pooled_cluster))]
  # adjusted_rand_indices <- sapply(dat %>% dplyr::select(ends_with("scSTMseq")), function(x) {
  #   adjustedRandIndex(x, dat$Group)
  # })
  
  # here we want to pull out the loglikelihood, both at 5, 10, 15, 30, 50 iteration, and the last iteration
  
  # ### SEURAT
  # seurat_name <- paste0(dirname(dirname(file_name)), "/seurat/seurat_",set_level,".rds")
  # seurat.sims <- safe_readRDS(seurat_name)
  # smeta <- seurat.sims@meta.data %>% as.data.frame()
  # sub_sims <- dat[match(rownames(smeta), dat$Cell),] # filter by the rows
  # seurat.adj <- sapply(smeta[,4:7], function(x) {
  #   adjustedRandIndex(x, sub_sims$Group)
  # })
  # # select the resolution that has the highest ARI. When there are multiple, select the first one
  # best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
  # seurat_ari <- max(seurat.adj)
  # names(seurat_ari) <- "seurat"

  # ##### calculating NMI ######
  # sce <- pooled_scSTMobj$settings$sce
  # x <- pooled_scSTMobj$settings$sce$Group
  # theta <- pooled_scSTMobj$theta
  # rownames(theta) <- pooled_scSTMobj$DocName
  # colnames(theta) <- paste0("Topic", 1:ncol(theta))
  # K <- pooled_scSTMobj$settings$dim$K
  # # calcluate entropy
  # collapsed.theta <- colSums(theta)
  # norm.collaposed.theta <- collapsed.theta/sum(collapsed.theta)
  # H.Y <- - sum(norm.collaposed.theta * log2(norm.collaposed.theta))
  # # calculate H(Y|C)
  # H.Y.C <- vector(mode = "list")
  # for(group in levels(x)){
  #   # extract cells that in the given group
  #   cells <- sce$Cell[sce$Group == group]
  #   subset.theta <- theta[rownames(theta) %in% cells,]
  #   collapsed.sub.theta <- colSums(subset.theta)
  #   norm.sub.theta <- collapsed.sub.theta/sum(collapsed.sub.theta)
  #   H.Y.C[[group]] <- - sum(norm.sub.theta * log2(norm.sub.theta))/K
  # }
  # H.Y.C <- do.call(sum, H.Y.C)
  # I.Y.C <- H.Y - H.Y.C
  # 
  # C.prop <- as.vector(prop.table(table(x)))
  # H.C <- - sum(C.prop*log2(C.prop))
  # NMI <- 2*I.Y.C/(H.Y + H.C)
  
  # Create a data frame to store results
  # res.temp <- data.frame(
  #   scSTMType = scSTMdir,
  #   # modelType = unlist(strsplit(set_level, "_"))[2],
  #   # seed = unlist(strsplit(set_level, "_"))[1],
  #   modelType = unlist(strsplit(set_level, "_"))[4],
  #   seed = unlist(strsplit(set_level, "_"))[3],
  #   nCellType = as.numeric(gsub("nCellType", "", unlist(strsplit(design, "_"))[2])),
  #   nSample = as.numeric(gsub("nSample", "", unlist(strsplit(design, "_"))[1])),
  #   Batch = ifelse(unlist(strsplit(design, "_"))[3]=="noBatch", FALSE, TRUE),
  #   CancerType = ifelse(unlist(strsplit(design, "_"))[4]=="StromalCell", FALSE, TRUE),
  #   Num_iter = as.numeric(gsub("iter","",unlist(strsplit(set_level, "_"))[1])),
  #   init_llh = pooled_scSTMobj$convergence$bound[1],
  #   Num_init = as.numeric(gsub("init","",unlist(strsplit(set_level, "_"))[2])),
  #   max_llh = tail(pooled_scSTMobj$convergence$bound, n=1),
  #   method = pooled_scSTMobj$settings$init$mode,
  #   llh5 = pooled_scSTMobj$convergence$bound[5],
  #   llh10 = pooled_scSTMobj$convergence$bound[10],
  #   llh15 = pooled_scSTMobj$convergence$bound[15],
  #   llh30 = pooled_scSTMobj$convergence$bound[30],
  #   llh50 = pooled_scSTMobj$convergence$bound[50],
  #   scSTMseq = t(adjusted_rand_indices)  # Transpose to match the original data frame structure
  #   # Seurat = seurat_ari
  # )
  # 
  # res <- bind_rows(res, res.temp)
  # cat(file_name, "\n")
  # rm(pooled_scSTMobj)
  # rm(dat)
}

write.csv(res, "res/1MultiSample_SingleResponse_Simulation/ARI_sensitivity_analysis_all.csv")


###################################################################################
############################ analysis #############################################
###################################################################################
res <- read.csv("res/1MultiSample_SingleResponse_Simulation/ARI_sensitivity_analysis_all.csv")
res$Num_iter <- factor(res$Num_iter, levels = as.numeric(names(table(res$Num_iter))))
res$Num_init <- factor(res$Num_init, levels = as.numeric(names(table(res$Num_init))))

# plot initualization vs final bound

seed_selection <- unique(res$seed)
ggplot(res %>% filter(seed == seed_selection[2] & Num_iter ==15), aes(x = init_llh, y = max_llh, color = method)) +
  geom_point() +
  facet_grid(~ seed, labeller = labeller(.cols = function(x) paste("Seed =", x))) +
  labs(title = "Bound at Initialization vs Final for a Single Simulation with Multiple Initialization",
       x = "Bound at Initialization",
       y = "Final Bound") +
  theme_bw()

ggplot(res %>% filter(seed == seed_selection[1] ), aes(x = init_llh, y = max_llh, color = Num_iter)) +
  geom_point() +
  facet_grid(~ method, labeller = labeller(.cols = function(x) paste("Method =", x)), scales = "free_x") +
  labs(title = "Bound at Initialization vs Final for a Single Simulation with Multiple Initialization",
       x = "Bound at Initialization",
       y = "Final Bound") +
  theme_bw()

ggplot(res  %>% filter(seed == seed_selection[1]) , aes(x = Num_iter, y = max_llh)) +
  geom_point() +
  facet_grid(~ Num_init, labeller = labeller(.cols = function(x) paste("Number of Initialization =", x))) +
  labs(x = "Number of Iterations", y = "ARI", title = "ARI vs Number of Iterations") +
  theme_bw()

###### plot final ARI versus final likelihood #########
ggplot(res, aes(x = pooled_scSTMseq, y = max_llh)) +
  geom_point() +
  facet_grid(~ Num_iter, labeller = labeller(.cols = function(x) paste("Number of Initial Intereation =", x))) +
  labs(x = "ARI", y = "Final Bound", title = "Final Bound versus ARI") +
  theme_bw()

###### likelihood at selection versus final likelihood ####

# pivot long to dynamicly plot llh
res_long <- res %>%
  pivot_longer(cols = starts_with("llh"), names_to = "loglikelihood", values_to = "llh_value")

res_long <- res_long %>%
  mutate(llh_num = as.numeric(gsub("llh", "", loglikelihood))) %>%
  filter(llh_num == Num_iter) # filter down to the llh that being selected

ggplot(res_long, aes(x = llh_value, y = max_llh)) +
  geom_point() +
  facet_grid(~ Num_iter, labeller = labeller(.cols = function(x) paste("Number of Initial Intereation =", x))) +
  labs(title = "Bound at Selection vs Final",
       x = "Bound at Selection",
       y = "Final Bound") +
  theme_minimal()


ggplot(res_long, aes(x = init_llh, y = max_llh, color = method)) +
  geom_point() +
  # facet_grid(~ Num_init, labeller = labeller(.cols = function(x) paste("Number of Initial States =", x))) +
  labs(title = "Bound at Initialization vs Final",
       x = "Bound at Initialization",
       y = "Final Bound") +
  theme_bw()


ggplot(res_long, aes(x = init_llh, y = max_llh, color = method)) +
  geom_point() +
  facet_grid(~ Num_init, labeller = labeller(.cols = function(x) paste("Number of Initial States =", x))) +
  labs(title = "Bound at Initialization vs Final",
       x = "Bound at Initialization",
       y = "Final Bound") +
  theme_bw()
cor(res$init_llh, res$max_llh)

ggplot(res_long %>% filter(Num_init == 5 ), aes(x = init_llh, y = max_llh)) +
  geom_point() +
  facet_grid(~ Num_init, labeller = labeller(.cols = function(x) paste("Number of Initial States =", x))) +
  labs(title = "Bound at Initialization vs Final",
       x = "Bound at Initialization",
       y = "Final Bound") +
  theme_bw()

res_long %>% filter(Num_init == 5 & init_llh < -8e8) 
# pre_res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final.csv")
# res <- read.csv("res/adjRandIndex_multiple_sample_benchmark_final_update.csv")
#
# res.update <- rbind(pre_res, res)
# res.update <- res %>%
#  # dplyr::select(sim, -matches("cluster$"), matches("cluster$")) %>%
#   dplyr::full_join(pre_res, by = c("sim", "seed", "control", "level", "nCellType", "nSample")) %>%
#  dplyr::select(sim, seed, control, level, nCellType, nSample,
#                -matches("cluster$"), matches("cluster$"))
# experiment_design <- sub(".*/SingleResponse/[^_]+_[^_]+_(.+)/.*", "\\1", paths[1])
# save_path <- paste0("res/1MultiSample_SingleResponse_Simulation/adjRandIndex_",basename(dir), "_", experiment_design, "_", basename(scSTMdir),".csv")
# write.csv(res, file = save_path)

# ##################################################################################
# ########################### plot #################################################
# ##################################################################################
# dat <- read.csv("res/1MultiSample_SingleResponse_Simulation/adjRandIndex_SingleResponse_Batch_CancerCell_scSTM_Pooled_Content_Sample_Prevalence_Time.csv")
#
# ari <- dat %>%
#   as.data.frame() %>%
#   dplyr::mutate(scSTMTypeNumeric = case_when(
#     modelType == "NullModel" ~ 0,
#     grepl("HighVar", modelType) ~ as.numeric(sub("HighVar", "", modelType))
#   ))
#
# # Create the box plot
# ggplot(ari, aes(x = factor(scSTMTypeNumeric), y = pooled_scSTMseq)) +
#   geom_boxplot() +
#   labs(
#     x = "scSTMType (0: NullModel, SD values)",
#     y = "ARI",
#     title = "ARI for Clustering Accuracy Using scSTMseq"
#   ) +
#   theme_minimal()
#
# dat_long <- dat %>%
#   mutate(nCellType = recode(nCellType, "5" = "5 Cell Types", "10" = "10 Cell Types",
#                             "15" = "15 Cell Types")) %>%
#   dplyr::rename("adjRandIndex" = "pooled_scSTMseq") %>%
#   gather(method, measurement, adjRandIndex:NMI, factor_key=TRUE) %>%
#   drop_na() %>%
#   group_by(nCellType, nSample, modelType, method) %>%
#   summarise(mean_measure = mean(measurement, na.rm = TRUE), .groups = "drop")
# # mutate(line_type = ifelse(grepl("scSTM", method), "solid", "dashed"))
#
# dat_long$nCellType <- factor(dat_long$nCellType, levels = c("5 Cell Types", "10 Cell Types", "15 Cell Types"))
# dat_long <- dat_long %>%
#   mutate(modelType = ifelse(grepl("Null",modelType), 0, as.numeric(sub("HighVar", "", modelType))))
# # dat_long$nSample <- factor(dat_long$nSample, levels = c("nsample3", "nsample6", "nsample12"))
#
# png("res/1MultiSample_SingleResponse_Simulation/Clustering_lineplot_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 1800, height = 1200, res = 220)
# ggplot(dat_long, aes(x = modelType, y = mean_measure, color = nCellType, group = nCellType)) +
#   # geom_line(aes(linetype = line_type), linewidth = 1.2) +
#   geom_line(linewidth = 1.2) +
#   labs(
#     x = "Effect Size",
#     y = "Mean Adjusted Rand Index",
#     title = "Mean Clustering Accuracy Across Cell Types and Effect Sizes"
#   ) +
#   facet_grid(~method) +
#   theme_bw() +
#   ylim(0, 1) +
#   theme(
#     text = element_text(size = 14),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     plot.title = element_text(size = 14, face = "bold"),
#     plot.subtitle = element_text(size = 12, face = "italic"),
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     panel.grid.major = element_line(color = "grey80"),
#     panel.grid.minor = element_blank(),
#     strip.text = element_text(size = 12, face = "bold")
#   )
# dev.off()
#
# dat_long <- dat %>%
#   mutate(nCellType = recode(nCellType, "5" = "5 Cell Types", "10" = "10 Cell Types",
#                             "15" = "15 Cell Types")) %>%
#   rename("adjRandIndex" = "pooled_scSTMseq") %>%
#   gather(method, measurement, adjRandIndex:NMI, factor_key=TRUE) %>%
#   drop_na()
# dat_long$nCellType <- factor(dat_long$nCellType, levels = c("5 Cell Types", "10 Cell Types", "15 Cell Types"))
# dat_long$modelType <- as.numeric(sub("HighVar", "", dat_long$modelType))
#
# png("res/1MultiSample_SingleResponse_Simulation/Clustering_boxplot_SingleResponse_noBatch_StromalCell_scSTM_Pooled_noContent_Prevalence_Time.png",
#     width = 1800, height = 1500, res = 220)
# ggplot(dat_long, aes(x = modelType, y = measurement, group = method)) +
#   geom_boxplot(outlier.shape = NA, width = 0.7) +
#   geom_jitter(size = 1, alpha = 0.6, aes(color = method)) +
#   facet_grid(method ~ nCellType) +
#   labs(
#     title = "Boxplot of Clustering Accuracy Across Cell Types and Effect Sizes",
#     x = "Effect Size",
#     y = "Measurements"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center and bold the title
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust x-axis text for better readability
#     axis.text.y = element_text(size = 12),  # Adjust y-axis text size
#     axis.title = element_text(size = 12),  # Increase axis title size
#     strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
#     legend.position = "none"  # Remove the legend as it's not necessary with the method labels on the x-axis
#   )
# dev.off()
