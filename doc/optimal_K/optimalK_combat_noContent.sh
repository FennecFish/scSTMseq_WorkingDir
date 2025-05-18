#!/bin/bash

#SBATCH --job-name=selectK
#SBATCH --output=selectK_%A_%a.out
#SBATCH --error=selectK_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-60
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/optimalK"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

# INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

# FILES=("sims_12_1714110786.rds" "sims_12_1714110802.rds" "sims_18_1714110807.rds" "sims_24_1714110767.rds" \
#  "sims_24_1714110854.rds" "sims_32_1714110781.rds" "sims_32_1714110810.rds" "sims_32_1714110827.rds" "sims_32_1714110843.rds" \
# "sims_5_1714110799.rds"  "sims_5_832594483.rds"   "sims_8_1714110755.rds"  "sims_8_1714110768.rds"  "sims_8_1714110849.rds")

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) 

Rscript optimalK_combat_noContent.R "${FILES[$INDEX]}"






