setwd("/proj/milovelab/wu/scLDAseq")
library(Matrix)
library(scran)
library(Rcpp)
library(slam)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stm)
library(scater)
library(scran)
library(doParallel)
library(foreach)
args <- commandArgs(trailingOnly = TRUE)
gc <- as.integer(args[1])
# dat <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/anti_PD1_cohort1_sce.rds")
# 
# # quick qc
# dat <- quickPerCellQC(dat, filter=TRUE)
# 
# ### remove genes with count 0 
# dat <- dat[rowSums(counts(dat)) != 0,]
# dat <- dat[,colSums(counts(dat)) != 0]
# #### feature selection #####
# dat <- scuttle::logNormCounts(dat)
# 
# dec.p2 <- modelGeneVar(dat)
# # feature selection
# p2.chosen <- getTopHVGs(dec.p2, n=2000)
# dat <- dat[p2.chosen,]
# cat("Feature Selection Completed")
# saveRDS(dat, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/PD1_cohort1_2000genes.rds")

dat <- readRDS("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_sce/PD1_cohort1_2000genes.rds")
dat$Batch <- paste0(dat$patient_id, "_", dat$expansion)
nsample <- length(unique(dat$patient_id))
ngroup <- length(unique(dat$cellType))

r.file <- paste0("R/",list.files("R/"))
sapply(r.file, source)
sourceCpp("src/STMCfuns.cpp")

t1 <- proc.time()
# scSTM.mod <- scSTMseq(sce = dat,
#                       sample = "patient_id", K = ngroup,
#                       prevalence = ~timepoint + patient_id, content = ~timepoint,
#                       init.type = "TopicScore",
#                       seed = 1, max.em.its = 20,
#                       verbose = TRUE)

# scSTM.mod <- scSTMseq(sce = dat, sample = "patient_id", K = ngroup,
#                                   prevalence = ~timepoint*expansion, content = ~Batch,
#                                   max.em.its = 5, init.type = "TopicScore", seed = 100)
# scSTM.mod <- selectModel_parallel(sce = dat, sample = "patient_id", K = ngroup,
#                           prevalence = ~timepoint*expansion, content = ~Batch,
#                           N = 5, ts_runs = 30, random_run =30,
#                          max.em.its = 100, net.max.em.its = 15, gc = gc)
scSTM.mod <- selectModel_parallel(sce = dat, sample = "patient_id", K = ngroup,
                                  prevalence = ~timepoint*expansion, content = ~Batch,
                                  N = 5, ts_runs = 30, random_run =30,
                                  interactions = FALSE,
                                  max.em.its = 100, net.max.em.its = 15, gc = gc)

# all_values <- unlist(scSTM.mod$bound)
# max_value <- max(all_values)
# max_position_in_vector <- which(all_values == max_value)
# res <- scSTM.mod$runout[[max_position_in_vector]]
msg <- sprintf("Completed scLDAseq (%d seconds). \n", floor((proc.time()-t1)[3]))
# saveRDS(scSTM.mod, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_Batch_Prevalence_TimeResponseInteraction.rds")
saveRDS(scSTM.mod, file = "/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_scSTM/scSTM_Content_BatchNoInteraction_Prevalence_TimeResponseInteraction.rds")