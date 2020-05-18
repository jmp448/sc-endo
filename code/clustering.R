rm(list = ls())
library(stats)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(tidyr)
library(superheat)
library(matrixStats)
library("viridis")           # Load

make.similarity <- function(my.data) {
  N <- nrow(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- cor(my.data[i,], my.data[j,], method = "pearson")
    }
  }
  S
}

make.affinity <- function(S, n.neighbors=2) {
  N <- length(S[,1])
  
  if (n.neighbors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighbors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}


metadata <- read.delim("/work-zfs/abattle4/prashanthi/sc-endo/data/metadata.pseudotime.tsv")
metadata$agg_group <- paste(metadata$donor_long_id, metadata$experiment, metadata$day, sep = "_")
pseudotime <- c()
for(isample in rownames(expr)){
  pseudotime <- c(pseudotime, mean(metadata$pseudo[metadata$agg_group == isample]))
}

expr <- readRDS("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/all_cells/expr.rds")
dynamic_eQTLs <- readRDS("/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/all_cells/stegle/significant_hits.rds")
dynamic_eQTLs$genes <- as.character(dynamic_eQTLs$genes)
dynamic_eQTLs$snps_interaction <- as.character(dynamic_eQTLs$snps_interaction)

expr <- expr[rownames(expr) %in% dynamic_eQTLs$genes, ]
metadata$agg_group <- paste(metadata$donor_long_id, metadata$experiment, metadata$day, sep = "_")
sample.id <- colnames(expr)
pseudotime <- c()
for(isample in colnames(expr)){
  pseudotime <- c(pseudotime, mean(metadata$pseudo[metadata$agg_group == isample]))
}

ordered.sample.id <- sample.id[order(pseudotime)]
expr <- expr[ ,match(ordered.sample.id, colnames(expr))]

S <- make.similarity(expr)

A <- make.affinity(S, 5)  # use 3 neighboors (includes self)
D <- diag(apply(A, 1, sum)) # sum rows
U <- D - A

"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))

k   <- 4
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]

set.seed(100)
km <- kmeans(Z, centers=k, nstart=5)
cluster.assignments <- km$cluster
genes <- rownames(expr)

ordered.genes <- genes[order(cluster.assignments)]
cluster.assignments <- cluster.assignments[order(cluster.assignments)]
genes_df <- cbind(ordered.genes, cluster.assignments)

expr <- expr[match(ordered.genes, rownames(expr)), ]

colSide <- brewer.pal(4, "Pastel1")[cluster.assignments[order(cluster.assignments)]]
sample.id <- colnames(expr)
day <- c()
for(id in sample.id){
  day <- c(day, strsplit(id, split = "_")[[1]][5])
}
pseudotime <- pseudotime[order(pseudotime)] 
pseudotime_factor <- pseudotime[order(pseudotime)]
pseudotime_factor[pseudotime < 0.15] <- 1
pseudotime_factor[(pseudotime >= 0.15) & (pseudotime <= 0.5)] <- 2
pseudotime_factor[(pseudotime > 0.5) & (pseudotime <= 0.7)] <- 3
pseudotime_factor[pseudotime > 0.7] <- 4

rowSide <- brewer.pal(4, "Set2")[as.factor(pseudotime_factor)]
pal <- colorRampPalette(brewer.pal(5, "RdBu"))(256)

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/expression_heatmap.pdf")
heatmap(-expr,  Rowv=NA, Colv=NA, col = pal, labRow = FALSE, labCol = FALSE,
        RowSideColors=colSide, ColSideColors = rowSide, margins = c(2, 2), main = "Expression clustering", xlab = "Pseudotime", ylab = "Genes")
dev.off()

cluster_1 <- colMeans(expr[cluster.assignments == 1, ])
cluster_2 <- colMeans(expr[cluster.assignments == 2, ])
cluster_3 <- colMeans(expr[cluster.assignments == 3, ])
cluster_4 <- colMeans(expr[cluster.assignments == 4, ])
df <- data.frame(pseudotime, cluster_1, cluster_2, cluster_3, cluster_4)

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/cluster_trajectories.pdf")
par(mfrow=c(2,2)) 
plot(cluster_1 ~ pseudotime, data = df, col = brewer.pal(4, "Pastel1")[1], pch = 20, ylab = "Expression", main = "Cluster 1")
with(df, lines(loess.smooth( pseudotime, cluster_1), col = "black", lwd = 0.8))

plot(cluster_2 ~ pseudotime, data = df, col = brewer.pal(4, "Pastel1")[2], pch = 20, ylab = "Expression", main = "Cluster 2")
with(df, lines(loess.smooth( pseudotime, cluster_2), col = "black", lwd = 0.8))

plot(cluster_3 ~ pseudotime, data = df, col = brewer.pal(4, "Pastel1")[3], pch = 20, ylab = "Expression", main = "Cluster 3")
with(df, lines(loess.smooth( pseudotime, cluster_3), col = "black", lwd = 0.8))

plot(cluster_4 ~ pseudotime, data = df, col = brewer.pal(4, "Pastel1")[4], pch = 20, ylab = "Expression", main = "Cluster 4")
with(df, lines(loess.smooth( pseudotime, cluster_4), col = "black", lwd = 0.8))
dev.off()

# effect size
effect <- matrix(rep(pseudotime,dim(dynamic_eQTLs)[1]) , nrow = dim(dynamic_eQTLs)[1], byrow = TRUE)
for(i in c(1:dim(dynamic_eQTLs)[1])){
  effect[i, ] <- effect[i, ]*dynamic_eQTLs$estimate[i]
}
cluster_1 <- dynamic_eQTLs$estimate[dynamic_eQTLs$genes %in% genes_df[genes_df[ ,2] == 1, 1]]
cluster_2 <- dynamic_eQTLs$estimate[dynamic_eQTLs$genes %in% genes_df[genes_df[ ,2] == 2, 1]]
cluster_3 <- dynamic_eQTLs$estimate[dynamic_eQTLs$genes %in% genes_df[genes_df[ ,2] == 3, 1]]
cluster_4 <- dynamic_eQTLs$estimate[dynamic_eQTLs$genes %in% genes_df[genes_df[ ,2] == 4, 1]]
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(256)
clusters <- c(rep(1, length(cluster_1)), rep(2, length(cluster_2)),
              rep(3, length(cluster_3)), rep(4, length(cluster_4)))
effect_df <- data.frame(clusters, c(cluster_1, cluster_2, cluster_3, cluster_4))
rownames(effect) <- dynamic_eQTLs$genes
colnames(effect) <- colnames(expr)
effect <- effect[match(ordered.genes, rownames(effect)), ]
pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/dynamic_model_predictions.pdf")
heatmap(effect, labRow = FALSE, labCol = FALSE, scale = "none", Colv = NA,
        ColSideColors = rowSide, margins = c(2, 2), 
        main = "eQTL effect", xlab = "Pseudotime", ylab = "Genes", col = coul)
dev.off()
# Read the deciles eQTLs
deciles <- list()
for(index in c(1:10)){
  coef_snps_by_gene <- list()
  if(index == 10){
    for(i in seq(0, 10584, 100)){
      start_pos <- i + 1
      end_pos <- min((start_pos + 99), 10584)
      coef_snps_by_gene <- c(coef_snps_by_gene, readRDS(paste0("../results/eQTL_calling/decile", as.character(index), "/10pc/coef_snps_by_gene_", start_pos, "to", end_pos, ".rds")))
  }}else{  
    for(i in seq(0, 10585, 100)){
    start_pos <- i + 1
    end_pos <- min((start_pos + 99), 10585)
    coef_snps_by_gene <- c(coef_snps_by_gene, readRDS(paste0("../results/eQTL_calling/decile", as.character(index), "/10pc/coef_snps_by_gene_", start_pos, "to", end_pos, ".rds")))
  }}

  nsnps <- lapply(coef_snps_by_gene, function(igene){
    igene <- as.data.frame(igene)
    dim(igene)[1]
  })
  nsnps <- unlist(nsnps)
  low.snps <- which(nsnps == 1)
  na_genes <- names(which(is.na(coef_snps_by_gene[low.snps])))
  coef_snps_by_gene <- coef_snps_by_gene[!names(coef_snps_by_gene) %in% na_genes]
  
  for(i in c(1:length(coef_snps_by_gene))){
    coef_snps_by_gene[[i]]$adj_p_value <- p.adjust(coef_snps_by_gene[[i]]$p.KR, method = "bonferroni")
  }
  coef_snps_by_gene <- coef_snps_by_gene[dynamic_eQTLs$genes]
  gene <- c()
  snps <- c()
  pvalue <- c()
  adj_pvalue <- c()
  for(j in c(1:dim(dynamic_eQTLs)[1])){
    gene[j] <- dynamic_eQTLs$genes[j]
    snps[j] <- dynamic_eQTLs$snps_interaction[j]
    a <- coef_snps_by_gene[[gene[j]]]
    if(!is.null(a)){
    pvalue[j] <- a$p.KR[rownames(a) == snps[j]]
    adj_pvalue[j] <- a$adj_p_value[rownames(a) == snps[j]]}
    else{
      pvalue[j] <- NA
      adj_pvalue[j] <- NA
    }
  }
  df <- data.frame(gene, snps, pvalue, adj_pvalue)
  colnames(df) <- c("gene", "snp", "pvalue", "adj_pvalue")
  deciles[[index]] <- df
}

dynamic_genes <- dynamic_eQTLs$genes
X <- matrix(data=NA,nrow=dim(dynamic_eQTLs)[1],ncol=10)
for(irow in c(1:length(dynamic_genes))){
  for(icol in c(1:10)){
    igene <- dynamic_genes[irow]
    a <- deciles[[icol]]
    if(igene %in% a$gene){
      X[irow, icol] <- -log10(a$adj_pvalue[a$gene == igene])}
    else{
      X[irow, icol] <- NA
    }
  }
}

rownames(X) <- dynamic_genes
colnames(X) <- paste0(rep("decile", 10), c(1:10))

X <- X[rowSums(is.na(X)) == 0, ]
X <- X[!rowSums(is.infinite(X)), ]
row_vars <- rowVars(X)
X <- X[!row_vars == 0, ]

genes_df <- genes_df[genes_df[ ,1] %in% rownames(X), ]
X <- X[match(genes_df[ ,1], rownames(X)), ]
pal.2 <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
colSide <- brewer.pal(4, "Pastel1")[as.factor(genes_df[ ,2])]

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/eQTL_clustering.pdf")
heatmap(-X,  Rowv=NA, Colv=NA, col = pal.2, labCol = FALSE, margins = c(2, 2), 
        main = "eQTL effect clustering", xlab = "Pseudotime", ylab = "Genes", RowSideColors=colSide,
        scale = "row", labRow = FALSE)
dev.off()

cluster_1_qtl <- colMeans(X[genes_df[ ,2] == 1, ])
cluster_2_qtl <- colMeans(X[genes_df[ ,2] == 2, ])
cluster_3_qtl <- colMeans(X[genes_df[ ,2] == 3, ])
cluster_4_qtl <- colMeans(X[genes_df[ ,2] == 4, ])
df_qtl <- data.frame(c(1:10), cluster_1_qtl, cluster_2_qtl, cluster_3_qtl, cluster_4_qtl)
colnames(df_qtl)[1] <- "decile"

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/eQTL_cluster_trajectories.pdf")
par(mfrow=c(2,2)) 
plot(cluster_1_qtl ~ c(1:10), data = df, col = brewer.pal(4, "Pastel1")[1], pch = 20, ylab = "-log10(p value)", main = "Cluster 1", xlab = "decile")
with(df, lines(loess.smooth( c(1:10), cluster_1_qtl), col = "black", lwd = 0.8))

plot(cluster_2_qtl ~ c(1:10), data = df, col = brewer.pal(4, "Pastel1")[2], pch = 20, ylab = "-log10(p value)", main = "Cluster 2", xlab = "decile")
with(df, lines(loess.smooth( c(1:10), cluster_2_qtl), col = "black", lwd = 0.8))

plot(cluster_3_qtl ~ c(1:10), data = df, col = brewer.pal(4, "Pastel1")[3], pch = 20, ylab = "-log10(p value)", main = "Cluster 3", xlab = "decile")
with(df, lines(loess.smooth( c(1:10), cluster_3_qtl), col = "black", lwd = 0.8))

plot(cluster_4_qtl ~ c(1:10), data = df, col = brewer.pal(4, "Pastel1")[4], pch = 20, ylab = "-log10(p value)", main = "Cluster 4", xlab = "decile")
with(df, lines(loess.smooth( c(1:10), cluster_4_qtl), col = "black", lwd = 0.8))
dev.off()
