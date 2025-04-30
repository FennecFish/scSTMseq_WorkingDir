#!/bin/bash

#SBATCH --job-name=adjRandIndex1000
#SBATCH --output=adjRandIndex1000.out
#SBATCH --error=adjRandIndex1000.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=30G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1
Rscript ARI_MultiResponse.R 

