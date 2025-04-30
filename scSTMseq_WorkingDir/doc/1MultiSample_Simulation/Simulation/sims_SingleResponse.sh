#!/bin/bash

#SBATCH --job-name=EffectSize_Sim
#SBATCH --output=EffectSize_Sim%A_%a.out
#SBATCH --error=EffectSize_Sim%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --array=1-20
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_SingleResponse.R $INDEX






