#!/bin/bash

#SBATCH --job-name=propeller
#SBATCH --output=propeller%A_%a.out
#SBATCH --error=propeller%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

Rscript propeller_replicates_eval.R 









