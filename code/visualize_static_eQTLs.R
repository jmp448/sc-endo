# Visualizing static eQTLs

rm(list = ls())
library(Biobase)
library(NMF)

res.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/"
deciles <- list()
all_hits <- list()
for(index in c(1:10)){
  deciles[[index]] <- readRDS(paste0(res.dir, "decile", as.character(index), "/10pc/significant_hits.rds"))
}

for(index in c(1:10)){
  all_hits[[index]] <- readRDS(paste0(res.dir, "decile", as.character(index), "/10pc/all_hits.rds"))
}


all.genes <- as.character(deciles[[1]]$genes)
for(index in c(2:10)){
  all.genes <- c(all.genes, as.character(deciles[[index]]$genes))
}

all.genes <- unique(all.genes)
X <- matrix(data=NA,nrow=length(all.genes),ncol=10)

for(irow in c(1:length(all.genes))){
  for(icol in c(1:10)){
    igene <- all.genes[irow]
    a <- all_hits[[icol]]
    if(igene %in% a$genes){
    X[irow, icol] <- -log10(a$adj_p_value[a$genes == igene])}
    else{
      X[irow, icol] <- NA
    }
  }
}

rownames(X) <- all.genes
colnames(X) <- paste0(rep("decile", 10), c(1:10))

X <- X[!rowSums(is.infinite(X)), ]

decile_correl <- cor(X)

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/decile_correlation.pdf")
superheat::superheat(decile_correl, pretty.order.rows = FALSE, 
                     pretty.order.cols = FALSE, left.label.col = "white",
                     bottom.label.col = "white", bottom.label.text.angle = 90,
                     bottom.label.text.size = 4, left.label.text.size = 4, 
                     grid.hline = FALSE, grid.vline = FALSE)
dev.off()

res <- nmf(X, 3)
pal <- colorRampPalette(brewer.pal(9, "Greens"))(256)

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/sparse_NMF.pdf")
heatmap(res@fit@H, Colv = NA, col = pal)
dev.off()

