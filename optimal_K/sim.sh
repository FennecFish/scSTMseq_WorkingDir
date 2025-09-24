#!/bin/bash

#SBATCH --job-name=sim_k
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-60
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

# Define the array of cell types
CellType=(5 8 12 18 24 32)

# Calculate the index for the job
index=$(($SLURM_ARRAY_TASK_ID / ${#CellType[@]}))

# Calculate the cell type index
cell_type_index=$(($SLURM_ARRAY_TASK_ID % ${#CellType[@]}))

# Get the cell type value from the array
CELLTYPE=${CellType[$cell_type_index]}

# Run the R script with the correct arguments
Rscript sim.R $index $CELLTYPE






