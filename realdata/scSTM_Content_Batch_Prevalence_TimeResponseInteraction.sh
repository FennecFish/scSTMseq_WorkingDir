#!/bin/bash

#SBATCH --job-name=PD1_Interaction
#SBATCH --output=PD1_Interaction.out
#SBATCH --error=PD1_Interaction.err
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --time=11-
#SBATCH --mem=800G
#SBATCH --array=1-2
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
for n in 4 6 8; do
  Rscript scSTM_Content_Batch_Prevalence_TimeResponseInteraction.R "$SLURM_NTASKS" $n
done








