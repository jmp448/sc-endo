#!/bin/bash -l
#SBATCH
#SBATCH --time=0:10:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=umap.err
#SBATCH --output=umap.out
#SBATCH --mail-type=END
#SBATCH --mail-user=jpopp4@jhu.edu

module load R/3.6.1

Rscript ../code/umap_significant.R

