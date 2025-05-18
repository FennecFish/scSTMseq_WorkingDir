#!/bin/bash

#SBATCH --job-name=PD1_effectSize
#SBATCH --output=PD1_effectSize.out
#SBATCH --error=PD1_effectSize.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1
Rscript MANOVA_eval.R
