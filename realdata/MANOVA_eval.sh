#!/bin/bash

#SBATCH --job-name=PD1_n4
#SBATCH --output=PD1_n4.out
#SBATCH --error=PD1_n4.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1
Rscript MANOVA_eval.R
