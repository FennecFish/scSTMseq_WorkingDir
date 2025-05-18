#!/bin/bash

#SBATCH --job-name=propellerBlockResponse
#SBATCH --output=propellerBlockResponse%A_%a.out
#SBATCH --error=propellerBlockResponse%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=5G
#SBATCH --array=11-962
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse/"
seurat_DIR="seurat"
FILES=($(find "$DIR" -type f -path "*/$seurat_DIR/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript propeller_block_sample_MultiResponse.R "$PARENT_DIR" "$seurat_DIR" "$FILE" 









