#!/bin/bash

#SBATCH --job-name=ManovaRM_nSample10_noBatch_Cancer
#SBATCH --output=ManovaRM_nSample10_noBatch_Cancer%A_%a.out
#SBATCH --error=ManovaRM_nSample10_noBatch_Cancer%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-4435
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
STM_DIR="scSTM_Pooled_noContent_Prevalence_Time"
FILES=($(find "$DIR" -type f -path "*nSample10_*_noBatch_Cancer*/$STM_DIR/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript ManovaRM_SingleResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







