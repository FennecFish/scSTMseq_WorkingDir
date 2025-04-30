#!/bin/bash

#SBATCH --job-name=Manova_pValue
#SBATCH --output=Manova_pValue%A_%a.out
#SBATCH --error=Manova_pValue%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=7G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1
Rscript pValue_Distribution_Under_Null.R





