#!/bin/bash

#SBATCH --job-name=scSTM_mixed
#SBATCH --output=scSTM_mixed%A_%a.out
#SBATCH --error=scSTM_mixed%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem=60G
#SBATCH --array=1-173
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

# ran 1-235
module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"

FILES=($(find "$DIR" -type f -path "*/sims/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript scSTM_Mixed_SingleResponse.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"







