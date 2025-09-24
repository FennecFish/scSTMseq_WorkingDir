#!/bin/bash

#SBATCH --job-name=scSTM_a_nc
#SBATCH --output=scSTM_a_nc_%A_%a.out
#SBATCH --error=scSTM_a_nc_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-
#SBATCH --mem-per-cpu=35G
#SBATCH --array=1-90
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_allGenes_noContent_run.R "${FILES[$INDEX]}"






