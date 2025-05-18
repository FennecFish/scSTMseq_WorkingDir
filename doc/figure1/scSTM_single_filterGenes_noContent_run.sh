#!/bin/bash

#SBATCH --job-name=scSTM_l3
#SBATCH --output=scSTM_l3_%A_%a.out
#SBATCH --error=scSTM_l3_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-202
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "*L3*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_single_filterGenes_noContent_run.R "${FILES[$INDEX]}"






