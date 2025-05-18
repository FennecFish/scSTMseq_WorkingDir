#!/bin/bash

#SBATCH --job-name=ttest
#SBATCH --output=ttest.out
#SBATCH --error=ttest.err
#SBATCH -n 1
#SBATCH --time=6:00:00
#SBATCH --mem=5g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript cc_single_ttest.R







