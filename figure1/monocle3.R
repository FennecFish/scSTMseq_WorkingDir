setwd("/proj/milovelab/wu/scLDAseq")
library(monocle3)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
file_name <- basename(file_name)
cat(file_name, "\n")

# file_name <- "sims_1712865827_L5.rds"
set_level <- sub("sims_([^.]*)\\.rds", "\\1",  file_name)


sims <- readRDS(paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims/", file_name))

expression_matrix <- counts(sims)
cell_metadata <- colData(sims)
gene_annotation <- rowData(sims)
gene_annotation$gene_short_name <- gene_annotation$Gene
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, num_dim = 100, alignment_group = "Batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
# k <- length(unique(sims$Group))
cds <- cluster_cells(cds, cluster_method = 'louvain')

saveRDS(cds, file = paste0("/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/monocle3/", 
                                   "monocle3_", set_level, ".rds"))
