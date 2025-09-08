#!/bin/bash

#SBATCH --job-name=ManovaRM_nSample10_noBatch_Cancer
#SBATCH --output=ManovaRM_nSample10_noBatch_Cancer%A_%a.out
#SBATCH --error=ManovaRM_nSample10_noBatch_Cancer%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-45
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

# need 302
module load r/4.3.1

SCSTM_DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/nSample10_nCellType5_noBatch_CancerCell/scSTM_Pooled_noContent_Prevalence_Time"
M_DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/nSample10_nCellType5_noBatch_CancerCell/ManovaRM"

# Find all .rds files in the scSTM directory
SCSTM_FILES=($(find "$SCSTM_DIR" -maxdepth 1 -type f -name "*.rds"))

# Initialize a list to hold scSTM files without matching ManovaRM counterparts
FILES=()

for FILE in "${SCSTM_FILES[@]}"; do
    # Extract base file name
    BASE_NAME=$(basename "$FILE")
    
    # Construct the expected ManovaRM file path by replacing the directory and filename prefix
    M_FILE="${M_DIR}/${BASE_NAME/scSTM_/ManovaRM_}"

    # If the corresponding ManovaRM file doesn't exist, save the scSTM file
    if [ ! -f "$M_FILE" ]; then
        FILES+=("$FILE")
    fi
done

# # Show all files in the FILES array
# for f in "${FILES[@]}"; do
#     echo "$f"
# done

# echo "There are ${#FILES[@]} files in total."

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

STM_DIR="scSTM_Pooled_noContent_Prevalence_Time"
echo $FILE
Rscript ManovaRM_SingleResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







