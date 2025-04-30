#!/bin/bash

#SBATCH --job-name=PD1_Interaction
#SBATCH --output=PD1_Interaction.out
#SBATCH --error=PD1_Interaction.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=4-
#SBATCH --mem=500G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript scSTM_Content_Batch_Prevalence_TimeResponseInteraction.R "$SLURM_NTASKS"







