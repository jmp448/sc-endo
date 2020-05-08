#!/bin/bash
#SBATCH --job-name=make_monocle
#SBATCH --time=1:0:0
#SBATCH --output=../log/make_monocle.out
#SBATCH --error=../log/make_monocle.err
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=bigmem2

module load R

Rscript ../code/make_monocle.R
