# Goal: This script is used to analyze scSTM from real data PD1
# 1) compare assigned cluster from paper cluster
# 2) changes between groups (time/expansion)
setwd("/proj/milovelab/wu/scLDAseq")
# setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(mclust)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

select_top_scSTM <- function(scSTMobj) {
  if(class(scSTMobj) == "selectModel") {
    all_values <- unlist(scSTMobj$bound)
    max_value <- max(all_values)
    if(length(which(all_values == max_value)) > 1){
      max_position_in_vector <- which(all_values == max_value)[1]
    } else{
      max_position_in_vector <- which(all_values == max_value)
    }
    scSTMobj <- scSTMobj$runout[[max_position_in_vector]]
  }
  return(scSTMobj)
}

cluster_scSTM <- function(scSTMobj) {
    max_indices <- apply(scSTMobj$theta, 1, which.max)
    colnames(scSTMobj$theta) <- paste0("topic_", 1:ncol(scSTMobj$theta))
    rownames(scSTMobj$theta) <- scSTMobj$DocName
    res_cluster <- colnames(scSTMobj$theta)[max_indices]
    names(res_cluster) <- rownames(scSTMobj$theta)
    return(res_cluster)
}


# scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")
scSTMobj <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_BatchNoInteraction_Prevalence_TimeResponseInteraction.rds")
scSTMobj <- select_top_scSTM(scSTMobj)
sims <- scSTMobj$settings$sce
# calculate adjusted rand index
scSTM_cluster <- cluster_scSTM(scSTMobj)
scSTM_cluster <- scSTM_cluster[match(rownames(colData(sims)), names(scSTM_cluster))]

adjustedRandIndex(scSTM_cluster, sims$cellType) # 0.6372073 for with batch interaction, 0.7622637 for without batch interaction

K <- length(unique(sims$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

png("res/PD1/structure_plot_Interaction.png", height = 2000, width = 4000, res = 300)
structure_plot(scSTMobj, topics = 1:K, grouping = sims$cellType, n = 3000, gap = 10)
dev.off()

# time effect
time_effect_i <- estimateEffect(1:K ~ timepoint, 
                              stmobj = scSTMobj, # ref = c("Pre", "NE"),
                              sampleNames = "patient_id",
                              sampleIDs = "BIOKEY_10",
                              uncertainty = "Global")
saveRDS(time_effect, file = "res/PD1/effect_Interaction.rds")

time_effect_ni <- estimateEffect(1:K ~ timepoint, 
                              stmobj = scSTMobj, # ref = c("Pre", "NE"),
                              sampleNames = "patient_id",
                              sampleIDs = "BIOKEY_10",
                              uncertainty = "Global")
# # try subset patients by expansion and non-expansion
# # expansion
# scSTMobj$settings$sce
# time_effect <- estimateEffect(1:K ~ timepoint, 
#                               stmobj = scSTMobj, # ref = c("Pre", "NE"),
#                               ssampleNames = "patient_id",
#                               sampleIDs = "BIOKEY_10",
#                               uncertainty = "Global")
# saveRDS(time_effect, file = "res/PD1/effect_Interaction.rds")
# summary(time_effect)

png("res/PD1/structure_plot_BatchNoInteraction.png", height = 800, width = 1200, res = 250)
structure_plot(scSTMobj, topics = 1:K, grouping = sims$cellType, n = 2000, gap = 20)
dev.off()

# time effect
ref <- c("Pre", "NE")
names(ref) <- c("timepoint", "expansion")
time_effect <- estimateEffect(1:K ~ timepoint*expansion,
                              stmobj = scSTMobj, ref.vec = ref,
                              uncertainty = "Global",
                              nsims = 30)
saveRDS(time_effect, file = "res/PD1/effect_BatchNoInteraction_nsim30.rds")

# 
# # time_effect <- readRDS("res/PD1/effect_BatchNoInteraction_nsim30.rds")
summary(time_effect)

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

png("res/PD1/estimateEffect_rep30.png", width = 1200, height = 800, res = 150)
plot(time_effect, covariate = "timepoint", model = scSTMobj,
     #method = "difference",
     method = "difference", # cov.value1="Pre",cov.value2="On",
     xlab = "proportion change between timepoints", moderator = "expansion",
     moderator.value = "E", linecol = "skyblue", printlegend = FALSE,
     main = "Expected Cell Group Proportion Change Before and After Treatment")

plot(time_effect, covariate = "timepoint", model = scSTMobj,
     #method = "difference",
     method = "difference", cov.value1="Pre",cov.value2="On",
     xlab = "timepoint", moderator = "expansion", 
     moderator.value = "NE", linecol = "red", add = T, 
     printlegend = T) 
# legend( c("expansion", "non-expansion"), lwd = 2, col = c("blue", "red"))



dat <- as.data.frame(t(labelTopics(scSTMobj, topics = 1:8, n = 50)$topics))
colnames(dat) <- paste0("Topic_", 1:8)
write.csv(dat, file = "res/PD1/top_50genes_in_topics.csv", row.names = F)

library(xtable)
latex_code <- xtable(dat)

labelTopics(scSTMobj, topics = 7, n = 5)$topics

###############################################################################
####################### Over-representation ###################################
###############################################################################
library(clusterProfiler)
library("AnnotationDbi")
library(org.Hs.eg.db)
n = 200
gene_int <- as.data.frame(t(labelTopics(scSTMobj, topics = 1:8, n = n)$topics))
colnames(gene_int) <- paste0("Topic_", 1:8)

gene_universe <- scSTMobj$vocab
# gene_universe <- labelTopics(scSTMobj, topics = 1:8, n = length(scSTMobj$vocab))
# gene_universe <- unique(as.vector(gene_universe$topics))

gene_universe <- unlist(
  mapIds(org.Hs.eg.db,
                     keys = gene_universe,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
  )

# total of 1915 
# remove 64 missing value
gene_universe <- gene_universe[!is.na(gene_universe)]

gene_lookup <- tibble(GeneSymbol = names(gene_universe), ENTREZID = gene_universe)

dat_long <- gene_int %>%
  pivot_longer(cols = everything(), names_to = "Topic", values_to = "GeneSymbol")

# Join with gene lookup table to replace gene symbols
dat_replaced <- dat_long %>%
  left_join(gene_lookup, by = "GeneSymbol") %>%
  # mutate(GeneSymbol = ifelse(is.na(ENTREZID), GeneSymbol, ENTREZID)) %>%
  dplyr::select(-GeneSymbol) %>%
  pivot_wider(names_from = Topic, values_from = ENTREZID) 
dat_replaced <- apply(dat_replaced, 2, unlist) 
dat_replaced <- split(dat_replaced, 1:ncol(dat_replaced))
dat_replaced <- lapply(dat_replaced, function(x){
  y <- x[!is.na(x)]
  return(y)
})
names(dat_replaced) <- paste0("Topic_", names(dat_replaced))
# compare clusters
ck <- compareCluster(geneCluster = dat_replaced, fun = "enrichKEGG")
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
png("res/PD1/KEGG_pathway.png", width = 4000, height = 4000, res = 400)
dotplot(ck, x="Cluster", showCategory = 10, title = "PD1 DataSet KEGG Pathway")
dev.off()

require(ReactomePA)
ck <- compareCluster(geneCluster = dat_replaced, fun = "enrichPathway")
png("res/PD1/reactome_pathway.png", width = 4000, height = 7000, res = 350)
dotplot(ck, x="Cluster", showCategory = 10, title = "PD1 DataSet Reactome Pathway")
dev.off()
# transfer ENTZ ID to Gene ID for plot
EzID <- as.character(ck@compareClusterResult$geneID)
geneID_readable <- sapply(EzID, function(x){
  ID <- gene_universe[which(gene_universe %in% unlist(strsplit(x, "/")))] #match each row back to Gene ID
  ID <- names(ID)
  return(paste0(ID,collapse = "/"))
})
geneID_readable <- unname(geneID_readable)
ck@compareClusterResult$geneID <- geneID_readable

png("res/PD1/reactome_pathway_network.png", width = 5000, height = 5000, res = 350)
cnetplot(ck)
dev.off()

ck <- compareCluster(geneCluster = dat_replaced, fun = "enrichDO")
png("res/PD1/enrichDO_pathway.png", width = 4000, height = 6000, res = 350)
dotplot(ck, x="Cluster", showCategory = 8)
dev.off()

ck <- compareCluster(geneCluster = dat_replaced, fun = "enrichGO", 
                     OrgDb='org.Hs.eg.db')
png("res/PD1/enrichGO_pathway.png", width = 4000, height = 7000, res = 300)
dotplot(ck, x="Cluster", showCategory = 10)
dev.off()
# disease ontology
# 
# ck <- compareCluster(geneCluster = dat_replaced, fun = "enrichWP", 
#                      organism = "human", pvalueCutoff = 0.1,
#                      qvalueCutoff = 1, pAdjustMethod = "none")


x <- enrichDO(gene          = dat_replaced$Topic_8,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = gene_lookup$ENTREZID,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

wiki <- enrichWP(as.numeric(dat_replaced$Topic_8), organism = "Homo sapiens") 
wiki@result[wiki@result$p.adjust < 0.1,]
# # ggo <- enrichGO(gene          = gene_list,
# #                 universe      = gene_universe,
# #                 OrgDb         = org.Hs.eg.db,
# #                 keyType       = "SYMBOL",
# #                 ont           = "BP",
# #                 pAdjustMethod = "BH",
# #                 pvalueCutoff  = 0.05,
# #                 qvalueCutoff  = 1,
# #                 readable      = TRUE)
# #
# # goplot(ggo)


#### KEGG
x <- enrichKEGG(
  gene = dat_replaced$Topic_7,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = gene_lookup$ENTREZID,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

x <- enrichPathway(
  gene = dat_replaced$Topic_7,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe =gene_lookup$ENTREZID,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
x@result[x@result$pvalue < 0.05,]

x <- lapply(dat_replaced, function(x){
  enrichWP(gene = x,
           organism = "Homo sapiens",
           universe =gene_lookup$ENTREZID,
           pvalueCutoff = 0.05,
           pAdjustMethod = "none",
           qvalueCutoff = 1)
})

dotplot(x$Topic_7, color = "pvalue")

x <- enrichWP(gene = dat_replaced$Topic_7,
              organism = "Homo sapiens",
              universe =gene_lookup$ENTREZID)
x@result[x@result$pvalue < 0.05,]


