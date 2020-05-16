rm(list = ls())

cell_type <- "all_cells/UMAP"
coef_snps_by_gene <- list()
pseudotime_by_gene <- list()
interaction_by_gene <- list()
for(i in seq(0, 10585, 100)){
  start_pos <- i + 1
  end_pos <- min((start_pos + 99), 10585)
  coef_snps_by_gene <- c(coef_snps_by_gene, readRDS(paste0("../results/eQTL_calling/", cell_type, "/10pc/coef_snps_by_gene_", start_pos, "to", end_pos, ".rds")))
  pseudotime_by_gene <- c(pseudotime_by_gene, readRDS(paste0("../results/eQTL_calling/", cell_type, "/10pc/pseudotime_by_gene_", start_pos, "to", end_pos, ".rds")))
  interaction_by_gene  <- c(interaction_by_gene , readRDS(paste0("../results/eQTL_calling/", cell_type, "/10pc/interaction_by_gene_", start_pos, "to", end_pos, ".rds")))
}
nsnps <- lapply(coef_snps_by_gene, function(igene){
  igene <- as.data.frame(igene)
  dim(igene)[1]
})
nsnps <- unlist(nsnps)
low.snps <- which(nsnps == 1)
na_genes <- names(which(is.na(coef_snps_by_gene[low.snps])))
coef_snps_by_gene <- coef_snps_by_gene[!names(coef_snps_by_gene) %in% na_genes]
pseudotime_by_gene <- pseudotime_by_gene[!names(pseudotime_by_gene) %in% na_genes]
interaction_by_gene <- interaction_by_gene[!names(interaction_by_gene) %in% na_genes]
for(i in c(1:length(coef_snps_by_gene))){
  coef_snps_by_gene[[i]]$adj_p_value <- p.adjust(coef_snps_by_gene[[i]]$p.KR, method = "bonferroni")
  pseudotime_by_gene[[i]]$adj_p_value <- p.adjust(pseudotime_by_gene[[i]]$p.KR, method = "bonferroni")
  interaction_by_gene[[i]]$adj_p_value <- p.adjust(interaction_by_gene[[i]]$p.KR, method = "bonferroni")
}
genes <- names(interaction_by_gene)
smallest_p_value_interaction <- c()
snps_interaction <- c()
genotype_pvalue <- c()
pseudotime_pvalue <- c()
for(igene in genes){
  isnp <- rownames(interaction_by_gene[[igene]])[which.min(interaction_by_gene[[igene]]$adj_p_value)]
  snps_interaction <- c(snps_interaction, isnp)
  smallest_p_value_interaction <- c(smallest_p_value_interaction, min(interaction_by_gene[[igene]]$adj_p_value))
  genotype_pvalue <- c(genotype_pvalue, coef_snps_by_gene[[igene]][isnp, ]$adj_p_value)
  pseudotime_pvalue <- c(pseudotime_pvalue, pseudotime_by_gene[[igene]][isnp, ]$adj_p_value)
}
final_df <- data.frame(genes, snps_interaction, smallest_p_value_interaction, genotype_pvalue, pseudotime_pvalue)
colnames(final_df) <- c("genes", "snps_interaction", "p_value_interaction", "p_value_genotype", "p_value_pseudotime")
final_df$adj_p_value_interaction <- p.adjust(final_df$p_value_interaction, method = "BH")
final_df$adj_p_value_genotype <- p.adjust(final_df$p_value_genotype, method = "BH")
final_df$adj_p_value_pseudotime <- p.adjust(final_df$p_value_pseudotime, method = "BH")
saveRDS(final_df, paste0("../results/eQTL_calling/", cell_type, "/10pc/all_hits.rds"))
sig_hits <- final_df[final_df$adj_p_value_interaction <= 0.1, ]
saveRDS(sig_hits, paste0("../results/eQTL_calling/", cell_type, "/10pc/significant_hits.rds"))
