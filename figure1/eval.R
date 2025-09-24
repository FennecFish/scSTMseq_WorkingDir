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
  max_indices <- apply(scSTMobj$theta, 1, which.max)
  colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
  rownames(scSTMobj$theta) <- colnames(scSTMobj$mu$mu)
  res_cluster <- colnames(scSTMobj$theta)[max_indices]
  names(res_cluster) <- rownames(scSTMobj$theta)
  return(res_cluster)
}

files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc/", pattern = "scSTM*")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/", pattern = "scSTM*")
dat <- data.frame()
lx <- paste0("L", 1:3)
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

    # filter gene with combat no content 
    combat_f_nc_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_f_nc/scSTM_",set_level,".rds")
    if(file.exists(combat_f_nc_name)){
      scSTM_combat_f_nc <- readRDS(combat_f_nc_name)
      scSTM_combat_f_nc_cluster <- process_scSTM(scSTM_combat_f_nc)
      dat$scSTM_combat_f_nc_cluster <- scSTM_combat_f_nc_cluster[match(dat$Cell, names(scSTM_combat_f_nc_cluster))]
    } else {dat$scSTM_combat_f_nc_cluster = NA}
    
    #STM
    stm_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/STM_f_nc/STM_",set_level,".rds")
    if(file.exists(stm_name)){
      stm.obj <- readRDS(stm_name)
      stm_cluster <- process_scSTM(stm.obj)
      dat$stm_cluster <- stm_cluster[match(dat$Cell, names(stm_cluster))]
    } else {dat$stm_cluster = NA}
    
    # Seurat
    seurat.sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/seurat/seurat_", set_level, ".rds"))
    smeta <- seurat.sims@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(smeta)] # filter by the rows
    seurat.adj <- sapply(smeta[,4:7], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(seurat.adj)[seurat.adj == max(seurat.adj)][1]
    seurat_cluster <- seurat.sims@meta.data %>% as.data.frame() %>% select(all_of(best_res))
    dat$seurat_cluster <- seurat_cluster[match(dat$Cell, rownames(seurat_cluster)),]
    rm(sub_sims)

    # sctransform
    sctf <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sctransform/sctransform_", set_level, ".rds"))
    sct_meta <- sctf@meta.data %>% as.data.frame()
    sub_sims <- sims[,rownames(sct_meta)] # filter by the rows
    sctf.adj <- sapply(sct_meta[,6:9], function(x) {
      adjustedRandIndex(x, sub_sims$Group)
    })
    # select the resolution that has the highest ARI. When there are multiple, select the first one
    best_res <- names(sctf.adj)[sctf.adj == max(sctf.adj)][1]
    sctf_cluster <- sctf@meta.data %>% as.data.frame() %>% select(all_of(best_res))
    dat$sctransform_cluster <- sctf_cluster[match(dat$Cell, rownames(sctf_cluster)),]
    rm(sub_sims)
    
    # fastTopics
    fasttopic_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/fastTopics/fastTopics_", set_level, ".rds")
    nmf.sims <- readRDS(fasttopic_name)
    max_indices <- apply(nmf.sims$L, 1, which.max)
    # colnames(nmf.sims$L) <- paste0("topic_", 1:ncol(nmf.sims$L))
    # rownames(nmf.sims$L) <- colnames(scSTMobj$mu$mu)
    fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    names(fastTopics_cluster) <- rownames(nmf.sims$L)
    dat$fastTopics_cluster <- fastTopics_cluster[match(dat$Cell, names(fastTopics_cluster))]
    # if(file.exists(fasttopic_name)){
    #   nmf.sims <- readRDS(fasttopic_name)
    #   max_indices <- apply(nmf.sims$L, 1, which.max)
    #   # colnames(nmf.sims$L) <- paste0("topic_", 1:ncol(nmf.sims$L))
    #   # rownames(nmf.sims$L) <- colnames(scSTMobj$mu$mu)
    #   fastTopics_cluster <- colnames(nmf.sims$L)[max_indices]
    #   names(fastTopics_cluster) <- rownames(nmf.sims$L)
    #   dat$fastTopics_cluster <- fastTopics_cluster[match(dat$Cell, names(fastTopics_cluster))]
    # } else {dat$fastTopics_cluster = NA}
    
    # monocle3
    monocle3 <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/monocle3/monocle3_", set_level, ".rds"))
    dat$monocle3_cluster <- partitions(monocle3)[match(dat$Cell, names(partitions(monocle3)))]
    
    # SC3
    sc3_file_name <- paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sc3/sc3_", set_level, ".rds")
    if(file.exists(sc3_file_name)){
      sc3 <- readRDS(sc3_file_name)
      sc3_cluster <- sc3@colData %>% as.data.frame() %>% select(ends_with("clusters"))
      dat$sc3_cluster <- sc3_cluster[match(dat$Cell, rownames(sc3_cluster)),]
    } else {dat$sc3_cluster = NA}
    
    adjusted_rand_indices <- sapply(dat[, 5:ncol(dat)], function(x) {
      adjustedRandIndex(x, sims$Group)
    })
    
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

write.csv(res.adj, file = "res/res_fig1_adj_sims_V3.csv")

# dat <- read.csv("res/res_fig1_adj_sims_V3.csv")
# dat <- dat %>% select(X,sim, level, scSTM_combat_f_nc_cluster) %>% filter(level %in% c("pos_L3", "neg_L3"))
# write.csv(dat, file = "res/scSTM_l3_adjRandIndex.csv")
# files <- list.files(path = "/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/scSTM_combat_f_nc/", 
#                     pattern = "scSTM_.*_pos_L3.*")
# 
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
