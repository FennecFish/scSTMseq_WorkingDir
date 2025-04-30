#!/bin/bash

#SBATCH --job-name=K40_mix-PD1
#SBATCH --output=K40_mix-PD1.out
#SBATCH --error=K40_mix-PD1.err
#SBATCH -n 1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript subset_anti-PD1_mix.R
#Rscript subset_PD1_mix.R


