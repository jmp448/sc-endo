#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=20:0:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/eQTL_calling_LMM.R $1 > ../log/$1.log
