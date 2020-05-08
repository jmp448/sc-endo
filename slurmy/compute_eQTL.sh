#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/calling_eQTL_LMM.R $1 $2 > ../log/$1_$2.log
