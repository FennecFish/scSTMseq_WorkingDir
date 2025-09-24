#!/bin/bash

#SBATCH --job-name=STM_Sensitivity_all
#SBATCH --output=STM_Sensitivity_all%A_%a.out
#SBATCH --error=STM_Sensitivity_all%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=8-
#SBATCH --mem=55G
#SBATCH --array=1-10
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
FILES=($(find "$DIR" -type f -path "*nSample10_nCellType8_noBatch*Cancer*/sims/*Null*.rds"| head -n 10))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

# Iterate over parameters
# for ITER in 5 30 50; do
#   for INIT in 5 15 50; do
# for ITER in 10 15; do
#   for INIT in 5 15 50; do
for ITER in 2; do
  for INIT in 30; do
    Rscript scSTM_Sensitivity.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS" "$ITER" "$INIT"
  done
done






