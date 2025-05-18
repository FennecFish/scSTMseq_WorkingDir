#!/bin/bash

#SBATCH --job-name=sim_MultiResponse
#SBATCH --output=sim_MultiResponse%A_%a.out
#SBATCH --error=sim_MultiResponse%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-10
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_MultiResponse.R $INDEX






