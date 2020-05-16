module load bedtools


for epigenome in `ls ../data/epigenomes/merged`
do
  for celltype in "all_cells/stegle" "all_cells/UMAP"
  do
    epi=`echo $epigenome | sed 's/.merged.bed//'`
    bedtools intersect -wb -a ../data/epigenomes/merged/$epigenome -b ../results/func_analysis/static.merged.bed  > ../results/func_analysis/$epi.static.bed
    echo done with $celltype/$epi.static.bed
    bedtools intersect -wb -a ../data/epigenomes/merged/$epigenome -b ../results/eQTL_calling/$celltype/10pc/significant_hits.bed  > ../results/func_analysis/$celltype/$epi.sig.bed
    echo done with ${celltype}/${epi}.bg.bed
  done
done
