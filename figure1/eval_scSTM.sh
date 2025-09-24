#!/bin/bash

#SBATCH --job-name=eval_scSTM
#SBATCH --output=eval_scSTM.out
#SBATCH --error=eval_scSTM.err
#SBATCH -n 1
#SBATCH --time=3:00:00
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval_scSTM.R






