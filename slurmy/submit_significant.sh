#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_sighits
#SBATCH --mail-type=END
#SBATCH --mail-user=jpopp4@jhu.edu
#SBATCH --time=0:5:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


for T in iPSC mesendo defendo day1 day3
do
    for C in 15
    do
	sbatch --job-name=sig.$T.$C --error=../log/errors/$T.$C.sig.err --output=../log/outputs/$T.$C.sig.out identify_significant.sh $T $C
    done
done
