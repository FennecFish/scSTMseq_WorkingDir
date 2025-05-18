#!/bin/bash

#SBATCH --job-name=sims_Null
#SBATCH --output=sims_Null%A_%a.out
#SBATCH --error=sims_Null%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --array=1-240
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_SingleResponse_Gamma.R $INDEX






