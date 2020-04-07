#!/bin/bash
#SBATCH --job-name=fig1
#SBATCH --time=30:0
#SBATCH --output=log/fig1.out
#SBATCH --error=log/fig1.err
#SBATCH --mem=40G
#SBATCH --partition=shared

module load python/3.7-anaconda

python fig1.py
