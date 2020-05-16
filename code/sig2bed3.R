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
get_chrom2 <- function(s) {
  c <- str_match(s, "(.*?)-")[2]
  return(c)
}
get_pos2 <- function(s) {
  c <- str_match(s, "--(.*)")[2]
  return(c)
}

find.closest <- function(a,bg) {
  i <- which.min(abs(a-bg))
  return(names(bg)[i])
}

#load genotype matrix
geno_matrix <- readRDS("../data/geno_matrix.rds")

# calculate MAF
maf <- rowSums(geno_matrix)/(2*ncol(geno_matrix))
maf <- sapply(maf, function(x){if(x>0.5){1-x} else {x}})

for (cell_type in c("all_cells/stegle", "all_cells/UMAP")) {
  
  # get snps
  sig <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/10pc/significant_hits.rds"))
  sig <- as.character(sig$snps)
  sig <- sig[sig!="geno_subset"]
  sig <- sig[!is.na(sig)]
  
  # write to bed file
  chrs.sig <- paste0("chr", sapply(sig, get_chrom))
  pos.sig <- as.numeric(sapply(sig, get_pos))
  starts.sig <- pos.sig-1
  bed.sig <- data.frame("chr"=chrs.sig, "start"=starts.sig, "end"=pos.sig)
  write.table(bed.sig, paste0("../results/eQTL_calling/", cell_type, "/10pc/significant_hits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
  # remove dashes
  sig <- str_remove_all(sig, "`")
  
  # get the mafs of the significant snps, and remove these from consideration
  maf.sig <- maf[sig]
  maf.bg <- maf[!names(maf) %in% sig]
  
  # take a random subset of 50000 snps
  rands <- sample(1:length(maf.bg), size=50000)
  maf.bg <- maf.bg[rands]
  
  # find the closest snps in maf.bg to the significant snps
  bg.snps <- rep(NA, length(maf.sig))
  for (i in 1:length(maf.sig)) {
    bg.snps[i] <- find.closest(maf.sig[i], bg=maf.bg)
    maf.bg <- maf.bg[names(maf.bg) != bg.snps[i]]
  }
  stopifnot(sum(is.na(bg.snps))==0)
  #bg.snps <- sapply(X=maf.sig, FUN=find.closest, bg=maf.bg)
  
  # save background snps to a bed file
  chrs.bg <- paste0("chr", sapply(bg.snps, get_chrom2))
  pos.bg <- as.numeric(sapply(bg.snps, get_pos2))
  starts.bg <- pos.bg-1
  bed.bg <- data.frame("chr"=chrs.bg, "start"=starts.bg, "end"=pos.bg)
  # sort bed file
  bed.bg <- bed.bg[order(bed.bg$chr, bed.bg$start),]
  write.table(bed.bg, paste0("../results/eQTL_calling/", cell_type, "/10pc/background_hits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

