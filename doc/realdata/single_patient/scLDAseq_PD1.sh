#!/bin/bash

#SBATCH --job-name=scLDAseq_PD1
#SBATCH --output=scLDAseq_PD1_%A_%a.out
#SBATCH --error=scLDAseq_PD1_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=20G
#SBATCH --array=2-30
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

CSV="../../patientID.csv"
PARAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CSV)
Rscript scLDAseq_PD1.R "$PARAM"




