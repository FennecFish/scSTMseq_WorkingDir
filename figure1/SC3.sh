#!/bin/bash

#SBATCH --job-name=sc3
#SBATCH --output=sc3_%A_%a.out.out
#SBATCH --error=sc3_%A_%a.out.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-612
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript SC3.R "${FILES[$INDEX]}"







