#!/bin/bash

#SBATCH --job-name=ManovaRM_nSample10_noBatch_Cancer
#SBATCH --output=ManovaRM_nSample10_noBatch_Cancer%A_%a.out
#SBATCH --error=ManovaRM_nSample10_noBatch_Cancer%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-455
#SBATCH --constraint=rhel8
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/"
STM_DIR="scSTM_Pooled_noContent_Prevalence_Time"
FILES=($(find "$DIR" -type f -path "*nSample10_nCellType5_noBatch_Cancer*/$STM_DIR/*.rds"))

# # Show all files in the FILES array
# for FILE in "${FILES[@]}"; do
#     echo "$FILE"
# done
# Display the count of files in the FILES array
# echo "There are ${#FILES[@]} files in total."


# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript ManovaRM_SingleResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







