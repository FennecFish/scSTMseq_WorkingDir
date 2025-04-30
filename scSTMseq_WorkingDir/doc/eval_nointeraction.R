# Goal: This script is used to analyze scSTM from real data PD1
# 1) compare assigned cluster from paper cluster
# 2) changes between groups (time/expansion)
# setwd("/proj/milovelab/wu/scLDAseq")
setwd("/Users/Euphy/Desktop/Research/Single_Cell_Cancer/scLDAseq")
library(Matrix)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(mclust)

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
scSTMobj_ni <- readRDS("data/scSTM_Content_TimeandResponse_Prevalence_TimeandSample.rds")
scSTMobj_ni <- select_top_scSTM(scSTMobj_ni)
sims_ni <- scSTMobj_ni$settings$sce
# calculate adjusted rand index
scSTM_ni_cluster <- cluster_scSTM(scSTMobj_ni)
scSTM_ni_cluster <- scSTM_ni_cluster[match(rownames(colData(sims_ni)), names(scSTM_ni_cluster))]
adjustedRandIndex(scSTM_ni_cluster, sims_ni$cellType)  # 0.615965

K <- length(unique(sims_ni$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)

png("res/PD1/structure_plot_noInteraction.png", height = 600, width = 1000, res = 200)
structure_plot(scSTMobj_ni, topics = 1:K, grouping = sims_ni$cellType, n = 2000, gap = 5)
dev.off()

# time effect
time_effect_ni <- estimateEffect(1:K ~ timepoint + expansion, 
                                stmobj = scSTMobj_ni, # ref = c("Pre", "NE"),
                                # sampleNames = "patient_id",
                                # sampleIDs = "BIOKEY_10",
                                uncertainty = "Global")
saveRDS(time_effect_ni, file = "res/PD1/effect_noInteraction.rds")

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
# plot.estimateEffect(time_effect, covariate = "timepoint", model=scSTMobj,
#                     method="difference",cov.value1="Pre",cov.value2="On", linecol = "black")

time_effect <- readRDS("res/PD1/effect_noInteraction.rds")
plot(time_effect, covariate = "timepoint", model = scSTMobj, 
     #method = "difference",
     method = "difference", cov.value1="Pre",cov.value2="On",
     xlab = "timepoint",linecol = "blue", printlegend = FALSE) 
# plot(time_effect, covariate = "timepoint", model = scSTMobj, 
#      #method = "difference",
#      method = "difference", cov.value1="Pre",cov.value2="On",
#      xlab = "timepoint", moderator = "expansion", 
#      moderator.value = "E", linecol = "blue", printlegend = FALSE) 

plot(time_effect, covariate = "timepoint", model = scSTMobj, 
     #method = "difference",
     method = "difference", cov.value1="Pre",cov.value2="On",
     xlab = "timepoint", moderator = "expansion", 
     moderator.value = "NE", linecol = "red", add = T, 
     printlegend = F) 
# legend( c("expansion", "non-expansion"), lwd = 2, col = c("blue", "red"))

labelTopics(scSTM_nC_P, 4)