#!/bin/bash

#SBATCH --job-name=scSTM_pooled
#SBATCH --output=scSTM_pooled%A_%a.out
#SBATCH --error=scSTM_pooled%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=3-
#SBATCH --mem=90G
#SBATCH --array=1-30
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

# for Multipresponse Sample 20
# ran 1-129, 500-550 (except 537,538)

module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/MultiResponse/"
FILES=($(find "$DIR" -type f -path "*/sims/*.rds"  -newermt "2024-11-26"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript scSTM_Pooled_MultiResponse.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"







