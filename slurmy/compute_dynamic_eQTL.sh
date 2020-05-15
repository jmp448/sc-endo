#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL
<<<<<<< HEAD
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:40:0
=======
#SBATCH --time=0:20:0
>>>>>>> beabdb18ad48d57e1fc7e9d53c465aed79849637
#SBATCH --partition=lrgmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load gcc/5.5.0
module load R/3.6.1

<<<<<<< HEAD
Rscript ../code/dynamic_eQTL_calling.R $1  > ../log/all_cells_quadratic/dynamic_$1.log
=======
Rscript ../code/dynamic_eQTL_calling.R $1 $2

>>>>>>> beabdb18ad48d57e1fc7e9d53c465aed79849637
