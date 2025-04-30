#!/bin/bash

#SBATCH --job-name=PD1_scSTM
#SBATCH --output=PD1_scSTM.out
#SBATCH --error=PD1_scSTM.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=5-
#SBATCH --mem=800G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript scSTM_Content_Patient_Prevalence_TimeandResponse.R "$SLURM_NTASKS"







