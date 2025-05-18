#!/bin/bash

#SBATCH --job-name=Manova
#SBATCH --output=Manova%A_%a.out
#SBATCH --error=Manova%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --array=1-18
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse/"
STM_DIR="scSTM_Pooled_Content_Sample_Prevalence_Time"
FILES=($(find "$DIR" -type f -path "*nSample40*/$STM_DIR/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript ManovaRM_MultiResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







