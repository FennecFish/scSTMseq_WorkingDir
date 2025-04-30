#!/bin/bash

#SBATCH --job-name=fastTopics
#SBATCH --output=fastTopics%A_%a.out
#SBATCH --error=fastTopics%A_%a.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=5G
#SBATCH --array=2-20
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
FILES=($(find "$DIR" -type f -path "*nSample10_nCellType10_Batch*Cancer*/sims/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript fastTopics.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"







