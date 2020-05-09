#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL
#SBATCH --mail-type=END
#SBATCH --mail-user=jpopp4@jhu.edu
#SBATCH --time=0:40:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/calling_eQTL_LMM.R $1 $2 $3 > ../log/$1_$2_$3.log
