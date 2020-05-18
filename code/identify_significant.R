rm(list = ls())
inputArgs <-  commandArgs(TRUE)
cell_type <- inputArgs[1]
npcs <- inputArgs[2]

coef_snps_by_gene <- list()
for(i in seq(0, 10584, 100)){
  start_pos <- i + 1
  end_pos <- min((start_pos + 99), 10584)
  coef_snps_by_gene <- c(coef_snps_by_gene, readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/coef_snps_by_gene_", start_pos, "to", end_pos, ".rds")))
}

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


genes <- names(coef_snps_by_gene)
smallest_p_value <- c()
snps <- c()
for(igene in genes){
  snps <- c(snps, rownames(coef_snps_by_gene[[igene]])[which.min(coef_snps_by_gene[[igene]]$adj_p_value)])
  smallest_p_value <- c(smallest_p_value, min(coef_snps_by_gene[[igene]]$adj_p_value))
}

final_df <- data.frame(genes, snps, smallest_p_value)
colnames(final_df) <- c("genes", "snps", "p_value")

final_df$adj_p_value <- p.adjust(final_df$p_value, method = "BH")

saveRDS(final_df, paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/all_hits.rds"))


significant_hits <- final_df[final_df$adj_p_value < 0.10, ]

saveRDS(significant_hits, paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/significant_hits.rds"))
