#!/bin/bash

#SBATCH --job-name=sim_single
#SBATCH --output=sim_single_%A_%a.out
#SBATCH --error=sim_single_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-7
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_single_sample_V3.R $INDEX





