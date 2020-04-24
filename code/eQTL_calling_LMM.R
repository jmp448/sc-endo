rm(list = ls())
library(nlme)
library(pbkrtest)

gene.locs <- read.delim("/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt", stringsAsFactors = FALSE)
gene.locs <- gene.locs[gene.locs$gene_type == "protein_coding", ]
gene.locs <- gene.locs[!gene.locs$chr == "chrM", ]
gene.locs <- gene.locs[!gene.locs$chr == "chrX", ]
gene.locs <- gene.locs[!gene.locs$chr == "chrY", ]
gene.locs$ensmbl <- sapply(strsplit(gene.locs$gene_id, split = ".", fixed = TRUE), `[`, 1)

cell_type <- "iPSC"

expr <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr.rds"))
expr <- t(expr)

genes <- data.frame(strsplit(colnames(expr), split = "_"))
colnames(genes) <- c(1:dim(genes)[2])
genes <- t(genes)
genes <- data.frame(genes)
genes$labels <- paste(genes$X1, genes$X2, sep = "_")
colnames(genes) <- c("ensmbl", "symbol", "labels")
genes$ensmbl <- as.character(genes$ensmbl)
gene.locs$ensmbl <- as.character(gene.locs$ensmbl)

genes <- genes[genes$ensmbl %in% gene.locs$ensmbl, ]
gene.locs <- gene.locs[gene.locs$ensmbl %in% genes$ensmbl, ]
gene.locs <- gene.locs[match(genes$ensmbl, gene.locs$ensmbl), ]
genes$chr <- gene.locs$chr
genes$start_pos <- gene.locs$start_pos
genes$end_pos <- gene.locs$end_pos
saveRDS(genes, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/gene_locs.rds"))

expr <- expr[ ,colnames(expr) %in% paste(genes$ensmbl, genes$symbol, sep = "_")]

# Read genotype and covariates
geno <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype.rds"))
geno <- t(geno)
sample_assignments <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/sample_assignments.rds"))
expr_PCs <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr_PCs.rds"))
geno_PCs <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype_PCs.rds"))

sample_assignments$sample_id <- as.character(sample_assignments$sample_id)

snps <- data.frame(strsplit(colnames(geno), split = "--"))
colnames(snps) <- c(1:dim(snps)[2])
snps <- t(snps)
snps <- data.frame(snps)
snps$labels <- colnames(geno)
colnames(snps) <- c("chr", "position", "labels")
snps$chr <- as.character(snps$chr)
snps$position <- as.character(snps$position)

snps$chr <- as.numeric(snps$chr)
snps$position <- as.numeric(snps$position)
rownames(snps) <- c(1:dim(snps)[1])
genes$chr <- as.numeric(gsub("chr","",genes$chr))

# 
window.size <- 250000
snps.select <- list()
for(i in c(1: dim(expr)[2])){
gene.id <- colnames(expr)[i]
chr <- genes$chr[genes$labels == gene.id]
start_pos <- genes$start_pos[genes$labels == gene.id]
end_pos <- genes$end_pos[genes$labels == gene.id]
start_pos <- start_pos - window.size
end_pos <- end_pos + window.size
isnps.select <- snps[snps$chr == chr, ]
isnps.select <- isnps.select[isnps.select$position >= start_pos, ]
isnps.select <- isnps.select[isnps.select$position <= end_pos, ]
snps.select[[i]] <- isnps.select
}

saveRDS(snps.select, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/snps_matched_genes.rds"))

colnames(expr_PCs) <- c("expr_PC1", "expr_PC2", "expr_PC3", "expr_PC4", 
                        "expr_PC5", "expr_PC6", "expr_PC7", "expr_PC8", 
                        "expr_PC9", "expr_PC10")

colnames(geno_PCs) <- c("geno_PC1", "geno_PC2", "geno_PC3", "geno_PC4", 
                        "geno_PC5", "geno_PC6", "geno_PC7")

# 
coef_snps_by_gene <- list()
for(i in c(1:dim(expr)[2])){
  print(i)
  geno_subset <- geno[ ,snps.select[[i]][ ,3]]
  print(dim(geno_subset))
  for(igeno in c(1:dim(geno_subset)[2])){
    fixed_eff <- cbind(geno_subset[ ,igeno], expr_PCs[ ,1:10], geno_PCs[ ,1:5])
    colnames(fixed_eff)[1] <- colnames(geno_subset)[igeno]
    data_df <- data.frame(expr[ ,i], fixed_eff, sample_assignments$individuals)
    colnames(data_df) <- c("expression", colnames(fixed_eff), "individuals")
    eQTL.model <- lmer(expression ~ (. - individuals) + (1 | individuals), data = data_df)
    icoefs <- data.frame(coef(summary(eQTL.model)))
    df.KR <- get_Lb_ddf(eQTL.model, fixef(eQTL.model))
    icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
    if(igeno == 1){
      coef_snps <- icoefs[2, ]
    }else{
      coef_snps <- rbind(coef_snps, icoefs[2, ])
    }
  }
  coef_snps_by_gene[[i]] <- coef_snps
}