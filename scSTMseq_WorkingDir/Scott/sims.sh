#!/bin/bash

#SBATCH --job-name=sim
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-10
#SBATCH --mail-type=all
#SBATCH --mail-user=youremail@live.unc.edu

# 10 arrays gives you 10 replicates
module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims.R $INDEX






