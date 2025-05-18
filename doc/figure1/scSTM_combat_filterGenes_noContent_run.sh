#!/bin/bash

#SBATCH --job-name=scSTM_combat_f_nc
#SBATCH --output=scSTM_combat_f_nc_%A_%a.out
#SBATCH --error=scSTM_combat_f_nc_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1376
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V1_multiple/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_combat_filterGenes_noContent_run.R "${FILES[$INDEX]}"






