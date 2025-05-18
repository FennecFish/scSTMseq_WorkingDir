#!/bin/bash

#SBATCH --job-name=STM_Null
#SBATCH --output=STM_Null%A_%a.out
#SBATCH --error=STM_Null%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=35G
#SBATCH --array=1-735
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
FILES=($(find "$DIR" -type f -path "*nSample10_*_noBatch*Cancer*/sims/*EffectSize*.rds" -newermt "2024-12-13"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript scSTM_Pooled_SingleResponse_noBatch.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"







