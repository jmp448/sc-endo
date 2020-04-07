#!/bin/bash
#SBATCH --job-name=fig1
#SBATCH --time=30:0
#SBATCH --output=slurm_output/fig1.out
#SBATCH --error=slurm_output/fig1.err
#SBATCH --mem=40G
#SBATCH --partition=shared
#SBATCH --mail-user=jpopp4@jhu.edu

module load python/3.7-anaconda

python fig1.py
