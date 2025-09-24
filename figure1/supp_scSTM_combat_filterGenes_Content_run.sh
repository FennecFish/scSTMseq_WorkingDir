#!/bin/bash

#SBATCH --job-name=L9_f_c
#SBATCH --output=L9_f_c_%A_%a.out
#SBATCH --error=L9_f_c_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-2
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

FILES=(
  sims_1715176272_neg_L1.rds
  sims_1715176272_pos_L1.rds
  sims_1715176286_neg_L1.rds
)


INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_combat_filterGenes_Content_run.R "${FILES[$INDEX]}"






