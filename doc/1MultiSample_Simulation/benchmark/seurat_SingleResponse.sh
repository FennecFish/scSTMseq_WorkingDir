#!/bin/bash

#SBATCH --job-name=seurat10_noBatch_Cancer
#SBATCH --output=seurat10_noBatch_Cancer%A_%a.out
#SBATCH --error=seurat10_noBatch_Cancer%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-4534
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
FILES=($(find "$DIR" -type f -path "*nSample10*_noBatch_CancerCell*/sims/*.rds"))
# FILES=($(find "$DIR" -type f -path "*/scSTM_Pooled_Content_Sample_Prevalence_Time/*.rds" -newermt "2024-12-04"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript seurat_SingleResponse.R "$PARENT_DIR" "$FILE"







