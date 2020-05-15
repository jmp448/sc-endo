rm(list = ls())
library("NMF")
res.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/"
iPSC <- readRDS(paste0(res.dir, "iPSC/significant_hits.rds"))
mesendo <- readRDS(paste0(res.dir, "mesendo/significant_hits.rds"))
defendo <- readRDS(paste0(res.dir, "defendo/significant_hits.rds"))
iPSC_mesendo <- intersect(iPSC$genes, mesendo$genes)
mesendo_defendo <- intersect(defendo$genes, mesendo$genes)
defendo_iPSC <- intersect(iPSC$genes, defendo$genes)
genes_common <- intersect(iPSC_mesendo, defendo$genes)
n.egenes <- c(dim(iPSC)[1], dim(mesendo)[1], dim(defendo)[1])
barplot(n.egenes,
main = "Number of eGenes (FDR 10%)",
xlab = "Cell population",
ylab = "eGenes",
names.arg = c("iPSC", "mesendo", "defendo"),
col = c("turquoise4", "tan1", "deeppink4"))
CPZ_defendo <- defendo[[grep("CPZ", names(defendo))]]
CPZ_iPSC <- iPSC[[grep("CPZ", names(iPSC))]]
CPZ_mesendo <- mesendo[[grep("CPZ", names(mesendo))]]
CPZ_defendo$snps <- rownames(CPZ_defendo)
CPZ_iPSC$snps <- rownames(CPZ_iPSC)
CPZ_mesendo$snps <- rownames(CPZ_mesendo)
CPZ_defendo$snps <- gsub("`4--", "", CPZ_defendo$snps)
CPZ_defendo$snps <- gsub("`", "", CPZ_defendo$snps)
CPZ_iPSC$snps <- gsub("`4--", "", CPZ_iPSC$snps)
CPZ_iPSC$snps <- gsub("`", "", CPZ_iPSC$snps)
CPZ_mesendo$snps <- gsub("`4--", "", CPZ_mesendo$snps)
CPZ_mesendo$snps <- gsub("`", "", CPZ_mesendo$snps)
CPZ_defendo$snps <- as.numeric(CPZ_defendo$snps)
CPZ_iPSC$snps <- as.numeric(CPZ_iPSC$snps)
CPZ_mesendo$snps <- as.numeric(CPZ_mesendo$snps)
plot(CPZ_iPSC$snps, -log10(CPZ_iPSC$p.KR), col = "turquoise4",
pch = 20, axes = F, xlab = "Base pair", ylab = "-log10(p value)", main = "CPZ: Lead Switching")
points(CPZ_mesendo$snps, -log10(CPZ_mesendo$p.KR), col = "tan1", pch = 20)
points(CPZ_defendo$snps, -log10(CPZ_defendo$p.KR), col = "deeppink4", pch = 20)
abline(v = 8621386, col = "turquoise4")
abline(v = 8567305, col = "deeppink4")
axis(1)
axis(2)
iPSC_common <- iPSC[iPSC$genes %in% genes_common,]
mesendo_common <- mesendo[mesendo$genes %in% genes_common, ]
defendo_common <- defendo[defendo$genes %in% genes_common, ]
common_hits <- data.frame(iPSC_common$adj_p_value, mesendo_common$adj_p_value, defendo_common$adj_p_value)
common_hits <- as.matrix(common_hits)
rownames(common_hits) <- iPSC_common$genes
common_hits <- -log10(common_hits)
common_hits <- common_hits[!is.infinite(rowSums(common_hits)), ]
heatmap(common_hits)
heatmap(common_hits, labRow = FALSE)
heatmap(common_hits, labRow = FALSE, cex = 0.5)
heatmap(common_hits, labRow = FALSE)
colnames(common_hits) <- c("iPSC", "mesendo", "defendo")
heatmap(common_hits, labRow = FALSE)
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(common_hits, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(common_hits, Colv = NA, Rowv = NA,col=colMain)
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(common_hits, labRow = FALSE, col = colMain)
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(common_hits, labRow = FALSE, col = colMain, cexCol = 0.6)
heatmap(common_hits, labRow = FALSE, col = colMain, cexCol = 0.8)
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(common_hits, labRow = FALSE, col = colMain, cexCol = 1.0)
