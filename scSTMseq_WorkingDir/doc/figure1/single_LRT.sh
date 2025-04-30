#!/bin/bash

#SBATCH --job-name=LRT
#SBATCH --output=LRT.out
#SBATCH --error=LRT.err
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=5g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript single_LRT.R







