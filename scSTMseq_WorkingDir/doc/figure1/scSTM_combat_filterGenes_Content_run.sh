#!/bin/bash

#SBATCH --job-name=scSTM_combat_f_c
#SBATCH --output=scSTM_combat_f_c_%A_%a.out
#SBATCH --error=scSTM_combat_f_c_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-
#SBATCH --mem-per-cpu=28G
#SBATCH --array=1-90
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_combat_filterGenes_Content_run.R "${FILES[$INDEX]}"






