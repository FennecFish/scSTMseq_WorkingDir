#!/bin/bash

#SBATCH --job-name=Manova_nsample10
#SBATCH --output=Manova_nsample10%A_%a.out
#SBATCH --error=Manova_nsample10%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=15G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1
Rscript eval_Manova_SingleResponse.R 






