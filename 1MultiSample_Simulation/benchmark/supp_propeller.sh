#!/bin/bash

#SBATCH --job-name=supp_propeller
#SBATCH --output=supp_propeller%A_%a.out
#SBATCH --error=supp_propeller%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=6G
#SBATCH --array=1-15
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

S_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse"

# Find all .rds files in the scSTM folder
S_FILES=($(find "$S_DIR" -type f -path "*/seurat/*.rds"))

# Filter out files that already have a corresponding propeller file starting with 'asin_' or 'logit_'
FILES=()
for SIM_FILE in "${S_FILES[@]}"; do
    # Extract the base name of the file
    BASE_NAME=$(basename "$SIM_FILE")
    
    # Identify the corresponding propeller directory path
    PROP_DIR=$(dirname "$SIM_FILE" | sed "s|seurat|propeller_unpaired|")

    # Construct the potential propeller file paths for both 'asin_' and 'logit_' prefixes
    ASIN_FILE="${PROP_DIR}/${BASE_NAME/seurat_/asin_}"
    LOGIT_FILE="${PROP_DIR}/${BASE_NAME/seurat_/logit_}"
    
    # Check if neither of the corresponding propeller files exists
    if [ ! -f "$ASIN_FILE" ] && [ ! -f "$LOGIT_FILE" ]; then
        FILES+=("$SIM_FILE")
    fi
done

# #Show all files in the FILES array
# for FILE in "${FILES[@]}"; do
#     echo "$FILE"
# done
# 
# # Display the count of files in the FILES array
# echo "There are ${#FILES[@]} files in total."

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")
STM_DIR="scSTM_Pooled_Content_Sample_Prevalence_Time"
seurat_DIR="seurat"
Rscript propeller_unpaired_MultiResponse.R "$PARENT_DIR" "$seurat_DIR" "$FILE" 


