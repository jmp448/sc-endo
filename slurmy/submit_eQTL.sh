#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_eQTL
#SBATCH --mail-type=END
#SBATCH --mail-user=jpopp4@jhu.edu
#SBATCH --time=0:5:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


for K in 0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000 5100 5200 5300 5400 5500 5600 5700 5800 5900 6000 6100 6200 6300 6400 6500 6600 6700 6800 6900 7000 7100 7200 7300 7400 7500 7600 7700 7800 7900 8000 8100 8200 8300 8400 8500 8600 8700 8800 8900 9000 9100 9200 9300 9400 9500 9600 9700 9800 9900 10000 10100 10200 10300 10400 10500
do
    for C in 10 15 20 30
    do
	if [ $K -eq 10500 ] && [ $C -eq 30 ]
	then
<<<<<<< HEAD
	    sbatch --job-name=eQTL_$1_$2_$3 --error=../log/errors/$1/$2-$3 --output=../log/outputs/$1/$2-$3 --mail-user=pravich2@jhu.edu compute_eQTL.sh $1 $K $C
=======
	    sbatch --job-name="eQTL-$1-$K-$C" --error="../log/errors/$1/$K-$C" --output="../log/outputs/$1/$K-$C" --mail-user=jpopp4@jhu.edu compute_eQTL.sh $1 $K $C
>>>>>>> beabdb18ad48d57e1fc7e9d53c465aed79849637
	else
	    sbatch --job-name="eQTL-$1-$K-$C" --error="../log/errors/$1/$K-$C" --output="../log/outputs/$1/$K-$C" compute_eQTL.sh $1 $K $C
	fi
    done
done
