module load bedtools
for epigenome in `ls ../data/epigenomes/merged`
do
  for celltype in "all_cells/stegle" "all_cells/UMAP"
  do
    epi=`echo $epigenome | sed 's/.merged.bed//'`
    bedtools intersect -wb -a ../data/epigenomes/merged/$epigenome -b ../results/eQTL_calling/$celltype/10pc/background_hits.bed  > ../results/func_analysis/$celltype/$epi.bg.bed
    echo done with $celltype/$epi.bg.bed
    bedtools intersect -wb -a ../data/epigenomes/merged/$epigenome -b ../results/eQTL_calling/$celltype/10pc/significant_hits.bed  > ../results/func_analysis/$celltype/$epi.sig.bed
    echo done with ${celltype}/${epi}.bg.bed
  done
done
