rm(list=ls())
library(stringr)

get_chrom <- function(s) {
  c <- str_match(s, "`(.*?)-")[2]
  return(c)
}
get_pos <- function(s) {
  c <- str_match(s, "--(.*?)`")[2]
  return(c)
}
for (cell_type in c("all_cells/stegle", "all_cells/UMAP")) {
  for (npcs in c(10)) {
	  #all <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/all_hits.rds"))
	  sig <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/significant_hits.rds"))
	  
	  #all <- as.character(all$snps)
	  #all <- all[all!="geno_subset"]
	  #all <- all[!is.na(all)]
	  sig <- as.character(sig$snps)
	  sig <- sig[sig!="geno_subset"]
	  sig <- sig[!is.na(sig)]
	  
	  #chrs.all <- paste0("chr", sapply(all, get_chrom))
	  #pos.all <- as.numeric(sapply(all, get_pos))
	  #starts.all <- pos.all-1
	  
	  chrs.sig <- paste0("chr", sapply(sig, get_chrom))
	  pos.sig <- as.numeric(sapply(sig, get_pos))
	  starts.sig <- pos.sig-1
	  
	  #bed.all <- data.frame("chr"=chrs.all, "start"=starts.all, "end"=pos.all)
	  bed.sig <- data.frame("chr"=chrs.sig, "start"=starts.sig, "end"=pos.sig)
	  
	  #write.table(bed.all, paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/all_hits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
	  write.table(bed.sig, paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/significant_hits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  }
}
