setwd("/work/users/e/u/euphyw/scLDAseq/data/PD1_data/PD1_raw")
library("Seurat")
library(SingleCellExperiment)
library(readr)
library(SeuratDisk)
set.seed(1)

################################################################################
############################## Data Process ####################################
################################################################################

# read in the count matrix and meta_data
cell1 <- readRDS("1863-counts_cells_cohort1.rds")
meta1 <- read_csv("1872-BIOKEY_metaData_cohort1_web.csv")
meta1 <- meta1[match(colnames(cell1),meta1$Cell),]

sce <- SingleCellExperiment(list(counts=cell1),
                            colData=DataFrame(meta1),
                            metadata=list(study="Anti_PD1"))
exclude_patient <- unique(colData(sce)$patient_id[colData(sce)$expansion=="n/a"])
sub_sce <- sce[, !sce$patient_id %in% exclude_patient]

saveRDS(sub_sce, "../PD1_sce/anti_PD1_cohort1_sce.rds")


cell2 <- readRDS("1867-counts_cells_cohort2.rds")
meta2 <- read_csv("1871-BIOKEY_metaData_cohort2_web.csv")

meta2 <- meta2[match(colnames(cell2),meta2$Cell),]

sce <- SingleCellExperiment(list(counts=cell2),
                            colData=DataFrame(meta2),
                            metadata=list(study="chemo"))
saveRDS(sub_sce, "../PD1_sce/chemo_cohort2_sce.rds")
