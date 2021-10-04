#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=core
#SBATCH --output=./outputs/%x-%j.out


# cd /home/ahalhed/AliciaMSc/squirrel
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome/scripts/core.R
