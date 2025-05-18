#!/bin/bash

#SBATCH --job-name=pd1_nobatch
#SBATCH --output=pd1_nobatch.out
#SBATCH --error=pd1_nobatch.err
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

Rscript sub_PD1_nobatch.R



