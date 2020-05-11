#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL.$1.$2.$3
#SBATCH --error=../log/errors/$1/$2-$3.err
#SBATCH --output=../log/outputs/$1/$2-$3.err
#SBATCH --time=1:0
#SBATCH --partition=express
if [ $2 == 10500 ] && [ $3 == 30 ] 
then
  #SBATCH --mail-user=jpopp4@jhu.edu
fi

module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/calling_eQTL_LMM.R $1 $2 $3 > ../log/log/$1/$2-$3.log

