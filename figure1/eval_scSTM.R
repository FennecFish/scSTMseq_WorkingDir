setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)
library(Seurat)
library(cluster)
library(monocle3)

process_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values)
    max_position_in_vector <- which(all_values == max_value)
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test5/", pattern = "scSTM*")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", pattern = "scSTM*")
dat <- data.frame()
lx <- paste0("L", 3)
control <- c("pos", "neg")
level <- expand.grid(lx, control) %>% as.data.frame() %>% mutate(level = paste0(Var2,"_",Var1))
level <- level$level
res.adj <- data.frame()

for(l in level){
  file_name <- grep(l, files, value = TRUE)
  # sim_name <- unique(sub(".*_([0-9]+_L[0-9]+).*", "\\1", file_name))
  res <- data.frame()
  
  for(file in file_name){
    set_level <- sub("scSTM_([^.]*)\\.rds", "\\1",  file)
    # set_level <- file
    sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims/sims_", set_level, ".rds"))
    dat <- colData(sims) %>% data.frame() %>% select(Cell:Group,time)
    
    scSTM_select <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test5/scSTM_",set_level,".rds")
    if(file.exists(scSTM_select)){
      scSTM_select <- readRDS(scSTM_select)
      scSTM_select_cluster <- process_scSTM(scSTM_select)
      dat$scSTM_select_cluster <- scSTM_select_cluster[match(dat$Cell, names(scSTM_select_cluster))]
    } else {dat$scSTM_select_cluster = NA}
    
    # # filter gene with combat no content 
    # scSTM_spectral <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_test3/scSTM_",set_level,".rds")
    # if(file.exists(scSTM_spectral)){
    #   scSTM_spectral <- readRDS(scSTM_spectral)
    #   scSTM_spectral_cluster <- process_scSTM(scSTM_spectral)
    #   dat$scSTM_spectral_cluster <- scSTM_spectral_cluster[match(dat$Cell, names(scSTM_spectral_cluster))]
    # } else {dat$scSTM_spectral_cluster = NA}
    # 
    # stm_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/STM_f_nc/STM_",set_level,".rds")
    # if(file.exists(stm_name)){
    #   stm.obj <- readRDS(stm_name)
    #   stm_cluster <- process_scSTM(stm.obj)
    #   dat$stm_cluster <- stm_cluster[match(dat$Cell, names(stm_cluster))]
    # } else {dat$stm_cluster = NA}
    
    # adjusted_rand_indices <- sapply(dat[, 5:ncol(dat)], function(x) {
    #   adjustedRandIndex(x, sims$Group)
    # })
    
    adjusted_rand_indices <- adjustedRandIndex(dat$Group, dat$scSTM_select_cluster)
    # Create a data frame to store results
    res.temp <- data.frame(
      sim = set_level,
      level = l,
      t(adjusted_rand_indices)  # Transpose to match the original data frame structure
    )
    # res.temp <- data.frame(
    #   sim = set_level,
    #   level = l,
    #   scSTM_a_nc_adjR = adjustedRandIndex(dat$scSTM_a_nc_cluster,sims$Group),
    #   scSTM_f_nc_adjR = adjustedRandIndex(dat$scSTM_f_nc_cluster,sims$Group),
    #   scSTM_f_c_adjR = adjustedRandIndex(dat$scSTM_f_c_cluster,sims$Group),
    #   Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
    #   fastTopics_adjR = adjustedRandIndex(dat$fastTopics_cluster,sims$Group))
    
    res <- bind_rows(res, res.temp)
    cat(file, "\n")
    rm(sims)
  }
  
  res.adj <- bind_rows(res.adj, res)
}

write.csv(res.adj, file = "res/scSTM_l3_adjRandIndex_t5.csv")
# 
# res <- read.csv("res/scSTM_l3_adjRandIndex_t5.csv")
# dat <- read.csv("res/scSTM_l3_adjRandIndex_t4.csv")
# dat$svd20 <- res$t.adjusted_rand_indices
# plot(dat$svd20, dat$scSTM_select_cluster)
# dat <- dat %>% select(X,sim, level, scSTM_combat_f_nc_cluster) %>% filter(level %in% c("pos_L3", "neg_L3"))
# write.csv(dat, file = "res/scSTM_l3_adjRandIndex.csv")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_combat_f_nc/",
#                     pattern = "scSTM_.*_pos_L3.*")

# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_combat_f_nc/", 
#                     pattern = "scSTM*")
# files <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_combat_f_nc/", files)
# 
# 
# res <- vector()
# for(file in files){
#   dat <- readRDS(file)
#   res[file] <- all(diff(dat$convergence$bound) > 0)
#   plot(dat$convergence$bound)
#   rm(dat)
# }
# read in scSTM
# 
# 
# msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# ###########################################################
# #################### Seurat ###############################
# ###########################################################
# t1 <- proc.time()
# 
# seurat.sims <- as.Seurat(sims, counts = "counts", data = "logcounts")
# seurat.sims <- FindVariableFeatures(seurat.sims, selection.method = "vst", nfeatures = 500)
# all.genes <- rownames(seurat.sims)
# seurat.sims <- ScaleData(seurat.sims, features = all.genes)
# seurat.sims <- RunPCA(seurat.sims)
# 
# seurat.sims <- FindNeighbors(seurat.sims, dims = 1:10)
# seurat.sims <- FindClusters(seurat.sims, resolution = 0.5)
# 
# 
# msg <- sprintf("Completed Seurat (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# saveRDS(seurat.sims, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/seurat_", sim_name, ".rds"))
# ###########################################################
# ###################### CIDR ###############################
# ###########################################################
# # t1 <- proc.time()
# # 
# # res.cidr <- scDataConstructor(as.matrix(counts(sims)))
# # res.cidr <- determineDropoutCandidates(res.cidr)
# # res.cidr <- wThreshold(res.cidr)
# # res.cidr <- scDissim(res.cidr)
# # res.cidr <- scPCA(res.cidr)
# # res.cidr <- nPC(res.cidr)
# # nCluster(res.cidr)
# # res.cidr <- scCluster(res.cidr)
# # dat$cidr_cluster <- res.cidr@clusters[match(colnames(res.cidr@tags), dat$Cell)]
# # 
# # msg <- sprintf("Completed CIDR (%d seconds). \n", floor((proc.time()-t1)[3]))
# # if(verbose) cat(msg)
# ###########################################################
# ##################### RACEID ###############################
# ###########################################################
# t1 <- proc.time()
# 
# # tutorial
# # https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
# sc <- SCseq(counts(sims))
# # Cells with a relatively low total number of transcripts are discarded.
# sc <- filterdata(sc,mintotal=2000) 
# # retrieve filtered and normalized expression matrix 
# # (normalized to the minimum total transcript count across all cells retained after filtering) 
# fdata <- getfdata(sc)
# # If all genes should be used, then the parameter FSelect needs to be set to FALSE. 
# sc <- compdist(sc, metric="pearson", FSelect = FALSE)
# # sc <- clustexp(sc)
# sc <- clustexp(sc,cln=ngroup,sat=FALSE) # FUNcluster for other methods
# sc <- findoutliers(sc)
# dat$raceID_cluster <- NA
# raceID_cluster <- sc@cluster$kpart[match(names(sc@cluster$kpart), dat$Cell)]
# dat$raceID_cluster <- sc@cluster$kpart[match(dat$Cell,names(sc@cluster$kpart))]
# msg <- sprintf("Completed RaceID (%d seconds). \n", floor((proc.time()-t1)[3]))
# if(verbose) cat(msg)
# 
# saveRDS(sc, file = paste0("/work/users/e/u/euphyw/scLDAseq/res/simulation/raceID_", sim_name, ".rds"))
# 
# set.seed(1)
# 
# #### evaluating all methods #####
# sc_eval <- function(sims, dat) {
#   
#   res <- data.frame()
#   # compute silhouette score
#   # dist.matrix <- dist(t(counts(sims)))
#   # scSTM.sil <- silhouette(as.numeric(as.factor(dat$scSTM_cluster)), dist.matrix)
#   # seurat.sil <- silhouette(as.numeric(as.factor(dat$seurat_cluster)), dist.matrix)
#   # raceid.sil <- silhouette(as.numeric(as.factor(dat$raceID_cluster)), dist.matrix)
#   # cidr.sil <- silhouette(as.numeric(as.factor(dat$cidr_cluster)), dist.matrix)
#   
#   res <- data.frame(
#     scSTM_adjR = adjustedRandIndex(dat$scSTM_cluster,sims$Group),
#     Seurat_adjR = adjustedRandIndex(dat$seurat_cluster, sims$Group), 
#     raceID_adjR = adjustedRandIndex(dat$raceID_cluster,sims$Group),#,
#     CIDR_adjR = adjustedRandIndex(dat$cidr_cluster,sims$Group))#,
#   # scSTM_sil = mean(scSTM.sil[,3]),
#   # seurat_sil = mean(seurat.sil[,3]),
#   # raceID_sil = mean(raceid.sil[,3]))
#   
#   return(res)
# }
# 
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/", pattern = "3cellTypes_sims.rds")
# sim_name <- sub("\\_sims.rds$", "", files)
# res <- data.frame()
# 
# for (i in sim_name){
#   dat <- read.csv(paste0("res/clustering_benchmark/colData_",i,".csv"))
#   sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",i,"_sims.rds"))
#   temp <- sc_eval(sims, dat = dat)
#   res <- rbind(res, temp)
#   rm(dat)
#   rm(sims)
#   cat(i, "\n")
# }
# 
# rownames(res) <- sim_name
# res <- tibble::rownames_to_column(res, "sim")
# res$NumCellType <- sub(".*samples_(\\d+)cellTypes", "\\1", res$sim)
# res$NumSample <- sub(".*_(\\d+samples).*", "\\1", res$sim)
# write.csv(res, file = "res/res_eval_clustering_4methods.csv")
# 
# # change the format from wide to long
# library(tidyr)
# res <- read.csv("res/res_eval_methods.csv")
# wcidr <- read.csv("res/res_eval_clustering_4methods.csv")
# 
# wcidr_long <- wcidr %>% gather(method, adjR, scSTM_adjR:CIDR_adjR, factor_key=TRUE)
# # gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# ggplot(wcidr_long, aes(x=method, y=adjR)) + 
#   geom_boxplot() + 
#   ggtitle("Comparison of Clustering Accuracy among Methods with 3 cell Types")
# 
# res <- res %>%
#   mutate(num_samples = as.numeric(gsub(".*_(\\d+)samples.*", "\\1", sim)),
#          num_cell_types = as.numeric(gsub(".*_(\\d+)cellTypes", "\\1", sim)))
# data_long <- res %>% filter(num_cell_types == 5) %>% gather(method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# # gather(res, method, adjR, scSTM_adjR:raceID_adjR, factor_key=TRUE)
# ggplot(data_long, aes(x=method, y=adjR)) + 
#   geom_boxplot() + 
#   ggtitle("Comparison of Clustering Accuracy among Methods with 5 cell Types")
# 
# 
# # looking at why scSTM is low
# 
# sim_name <- "BIOKEY_11_10samples_5cellTypes"
# sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/",sim_name,"_sims.rds"))
