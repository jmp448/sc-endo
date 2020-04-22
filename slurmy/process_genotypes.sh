#!/bin/bash -l

homeDir=/work-zfs/abattle4/prashanthi/sc-endo
rawDir=$homeDir/data/genotypes/raw_genotyping_array_calls


module load bcftools
module load vcftools
module load htslib

rawfiles=`ls $rawDir/*.vcf`
for file in $rawfiles
do
	baseFilename=`basename $file *.vcf`
	echo $baseFilename
	bgzip -c $file > $file.gz
	tabix -p vcf $file.gz
	for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
	do
		processedDir=$homeDir/data/genotypes/chr$chr
		tabix -h $file.gz $chr > $processedDir/$baseFilename
	done
done

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	processedDir=$homeDir/data/genotypes/chr$chr
	for file in $processedDir/*
	do
		bgzip -c $file > $file.gz
		tabix $file.gz
	done
	bcftools merge $processedDir/HPS*.vcf.gz -Oz -o $processedDir/all_individuals.vcf.gz
	vcftools --gzvcf $processedDir/all_individuals.vcf.gz --maf 0.05 --max-maf 0.5 --hwe 1e-8 --012 --max-missing 0.95 --out $processedDir/combined.chr$chr.common
done

