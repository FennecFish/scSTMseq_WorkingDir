#!/bin/bash

#SBATCH --job-name=fastTopics
#SBATCH --output=fastTopics_%A_%a.out
#SBATCH --error=fastTopics_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-612
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript fastTopics.R "${FILES[$INDEX]}"






