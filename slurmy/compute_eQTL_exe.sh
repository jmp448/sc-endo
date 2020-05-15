
module load gcc/5.5.0
module load R/3.6.1

Rscript ../code/calling_eQTL_LMM.R $1 $2 $3 > ../log/$1/$3pc/start_$2.log

