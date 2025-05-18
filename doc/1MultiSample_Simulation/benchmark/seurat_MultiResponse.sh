#!/bin/bash

#SBATCH --job-name=seurat
#SBATCH --output=seurat_%A_%a.out
#SBATCH --error=seurat_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --array=151-1008
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse/"
FILES=($(find "$DIR" -type f -path "*/sims/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

module load r/4.3.1
Rscript seurat_MultiResponse.R "$PARENT_DIR" "$FILE"







