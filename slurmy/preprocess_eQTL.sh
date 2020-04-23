#!/bin/bash -l
#SBATCH
#SBATCH --job-name=glasso
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=10:0:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/eQTL_calling_preprocessing.R $1
