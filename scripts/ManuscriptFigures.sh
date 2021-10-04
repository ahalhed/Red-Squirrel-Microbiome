#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:15:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=Figures
#SBATCH --output=./outputs/%x-%j.out


module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome/scripts/ManuscriptFigures.R
