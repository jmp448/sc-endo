#!/bin/bash -l
#SBATCH
#SBATCH --time=0:10:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load R/3.6.1

Rscript ../code/identify_significant.R $1 $2 > ../log/log/$1/$2-sighits.log

