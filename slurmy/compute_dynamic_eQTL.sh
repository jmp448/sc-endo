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

Rscript ../code/dynamic_eQTL_calling.R $1  > ../log/all_cells_quadratic/dynamic_$1.log
