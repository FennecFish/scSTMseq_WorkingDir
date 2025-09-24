#!/bin/bash

#SBATCH --job-name=sp_scSTM_a_c
#SBATCH --output=sp_scSTM_a_c_%A_%a.out
#SBATCH --error=sp_scSTM_a_c_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-
#SBATCH --mem-per-cpu=150G
#SBATCH --array=1-12
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

FILES=(
sims_1712873779_L4.rds
sims_1712873779_L5.rds
sims_1712873779_L7.rds
sims_1712873779_L8.rds
sims_1712873781_L3.rds
sims_1712873781_L4.rds
sims_1712873781_L5.rds
sims_1712873781_L6.rds
sims_1712873781_L8.rds
sims_1712873781_L9.rds
sims_1712873792_L8.rds
sims_1712873792_L9.rds
)

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_allGenes_Content_run.R "${FILES[$INDEX]}"






