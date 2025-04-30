#!/bin/bash

#SBATCH --job-name=eval_methods
#SBATCH --output=eval_methods.out
#SBATCH --error=eval_methods.err
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R






