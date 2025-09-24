#!/bin/bash

#SBATCH --job-name=sp_scSTM_f_c
#SBATCH --output=sp_scSTM_f_c_%A_%a.out
#SBATCH --error=sp_scSTM_f_c_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=100G
#SBATCH --array=1-5
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

FILES=(
    sims_1712873833_L7.rds
    sims_1712873832_L9.rds
    sims_1712873832_L8.rds
    sims_1712873832_L2.rds
    sims_1712873811_L9.rds
)


INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_filterGenes_Content_run.R "${FILES[$INDEX]}"






