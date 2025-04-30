#!/bin/bash

#SBATCH --job-name=stm_PD1
#SBATCH --output=stm_PD1_%A_%a.out
#SBATCH --error=stm_PD1_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem-per-cpu=40G
#SBATCH --array=2-30
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

CSV="../../patientID.csv"
PARAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CSV)
Rscript stm_anti-PD1.R "$PARAM"




