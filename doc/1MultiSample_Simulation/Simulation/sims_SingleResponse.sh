#!/bin/bash

#SBATCH --job-name=CellType5
#SBATCH --output=CellType5%A_%a.out
#SBATCH --error=CellType5%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-100
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_SingleResponse.R $INDEX






