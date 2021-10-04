#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=SU08-month
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "05/22/2020"

#script starts here
#---

#set up
# cd /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace AG08 with specific grid/year combo being run
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/scripts/month-SU08.R
