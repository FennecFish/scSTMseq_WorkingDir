#!/bin/bash

#SBATCH --job-name=scSTM_nCell5
#SBATCH --output=scSTM_nCell5%A_%a.out
#SBATCH --error=scSTM_nCell5%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=38G
#SBATCH --constraint=rhel8
#SBATCH --array=1-322
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.4.0

# DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/"
# FILES=($(find "$DIR" -type f -path "*nSample10_nCellType5_noBatch_CancerCell*/sims/*.rds"))
SIMS_DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/nSample10_nCellType5_noBatch_CancerCell/sims"
SCSTM_DIR="/users/e/u/euphyw/scSTMseq/data/simulation/1MultiSample/SingleResponse/nSample10_nCellType5_noBatch_CancerCell/scSTM_Pooled_noContent_Prevalence_Time"

SIMS_FILES=($(find "$SIMS_DIR" -maxdepth 1 -type f -name "*"))
# Filter out files that already have a corresponding 'scSTM_' file
FILES=()
for SIM_FILE in $(find $SIMS_DIR -type f -name "*.rds"); do
    # Extract the base name of the file
    BASE_NAME=$(basename "$SIM_FILE")
    
    # Identify the corresponding scSTM directory path by substituting 'sims' with 'scSTM_Pooled_noContent_Prevalence_Time'
    DIRECTORY=$(dirname "$SIM_FILE" | sed "s|sims|scSTM_Pooled_noContent_Prevalence_Time|")

    # Construct the corresponding scSTM file name
    SCSTM_FILE="${DIRECTORY}/${BASE_NAME/sims_/scSTM_}"
    
    # Check if the corresponding scSTM file does not exist
    if [ ! -f "$SCSTM_FILE" ]; then
        FILES+=("$SIM_FILE")
    fi
done

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

Rscript scSTM_Pooled_SingleResponse_noBatch.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"







