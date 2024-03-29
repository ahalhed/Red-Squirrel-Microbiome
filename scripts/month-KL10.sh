#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=KL10-month
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "10/04/2021"

#script starts here
#---

#set up
# cd /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome/scripts/month-KL10.R
