#!/bin/bash

#SBATCH --job-name=PD1_TimeandResponse
#SBATCH --output=PD1_TimeandResponse.out
#SBATCH --error=PD1_TimeandResponse.err
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --time=11-
#SBATCH --mem=700G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript scSTM_Content_TimeandResponse_Prevalence_TimeandSample.R "$SLURM_NTASKS"







