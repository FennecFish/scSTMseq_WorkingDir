#!/bin/bash

#SBATCH --job-name=Manova1000
#SBATCH --output=Manova%A_%a.out
#SBATCH --error=Manova%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-232
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
STM_DIR="scSTM_Pooled_Content_Sample_Prevalence_Time"
FILES=($(find "$DIR" -type f -path "*_Batch_*/$STM_DIR/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript ManovaRM_Max_SingleResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







