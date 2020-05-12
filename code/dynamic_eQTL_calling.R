rm(list = ls())
library(nlme)
library(pbkrtest)

cell_type <- "all_cells"
inputArgs <-  commandArgs(TRUE)

expr <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr.rds"))
expr <- t(expr)

genes <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/gene_locs.rds"))
expr <- expr[ ,colnames(expr) %in% paste(genes$ensmbl, genes$symbol, sep = "_")]

geno <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype.rds"))
geno <- t(geno)
sample_assignments <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/sample_assignments.rds"))
expr_PCs <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr_PCs.rds"))
geno_PCs <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype_PCs.rds"))

sample_assignments$sample_id <- as.character(sample_assignments$sample_id)

snps.select <- readRDS(paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/snps_matched_genes.rds"))

colnames(expr_PCs) <- c("expr_PC1", "expr_PC2", "expr_PC3", "expr_PC4", 
                        "expr_PC5", "expr_PC6", "expr_PC7", "expr_PC8", 
                        "expr_PC9", "expr_PC10", "expr_PC11", "expr_PC12", 
                        "expr_PC13", "expr_PC14", "expr_PC15")

colnames(geno_PCs) <- c("geno_PC1", "geno_PC2", "geno_PC3", "geno_PC4", 
                        "geno_PC5", "geno_PC6", "geno_PC7")

start_pos <- as.numeric(inputArgs[1]) + 1
#end_pos <- dim(expr)[2]
end_pos <- min(start_pos + 99, dim(expr)[2])

# Read in the pseudotime 
metadata <- read.delim("/work-zfs/abattle4/prashanthi/sc-endo/data/metadata.pseudotime.tsv")
metadata$agg_group <- paste(metadata$donor_long_id, metadata$experiment, metadata$day, sep = "_")
pseudotime <- c()
for(isample in rownames(expr)){
  pseudotime <- c(pseudotime, mean(metadata$pseudo[metadata$agg_group == isample]))
}

coef_snps_by_gene <- list()
pseudotime_by_gene <- list()
interaction_by_gene <- list()
index <- 1
expr <- scale(expr)
for(i in c(start_pos:end_pos)){
  print(i)
  geno_subset <- geno[ ,snps.select[[i]][ ,3]]
  geno_subset <- as.data.frame(geno_subset)
  print(dim(geno_subset))
  if(dim(geno_subset)[2] > 0){
    for(igeno in c(1:dim(geno_subset)[2])){
      fixed_eff <- cbind(geno_subset[ ,igeno], pseudotime, pseudotime*geno_subset[ ,igeno], expr_PCs[ ,1:10], geno_PCs[ ,1:5])
      colnames(fixed_eff)[1] <- colnames(geno_subset)[igeno]
      colnames(fixed_eff)[3] <- "Gxt"
      data_df <- data.frame(expr[ ,i], fixed_eff, sample_assignments$individuals)
      colnames(data_df) <- c("expression", colnames(fixed_eff), "individuals")
      eQTL.model <- lmer(expression ~ (. - individuals) + (1 | individuals), data = data_df)
      icoefs <- data.frame(coef(summary(eQTL.model)))
      df.KR <- get_Lb_ddf(eQTL.model, fixef(eQTL.model))
      icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
      if(igeno == 1){
        coef_snps <- icoefs[2, ]
        coef_time <- icoefs[3, ]
        coef_interaction <- icoefs[4, ]
      }else{
        coef_snps <- rbind(coef_snps, icoefs[2, ])
        coef_time <- rbind(coef_time, icoefs[3, ])
        coef_interaction <- rbind(coef_interaction, icoefs[4, ])
      }
      rownames(coef_time) <- rownames(coef_snps)
      rownames(coef_interaction) <- rownames(coef_snps)
    }
    coef_snps_by_gene[[index]] <- coef_snps
    pseudotime_by_gene[[index]] <- coef_time
    interaction_by_gene[[index]] <- coef_interaction
  }else{
    coef_snps_by_gene[[index]] <- NA
    pseudotime_by_gene[[index]] <- NA
    interaction_by_gene[[index]] <- NA
  }
  index <- index + 1
}

names(coef_snps_by_gene) <- names(snps.select)[start_pos:end_pos]
names(pseudotime_by_gene) <- names(snps.select)[start_pos:end_pos]
names(interaction_by_gene) <- names(snps.select)[start_pos:end_pos]

saveRDS(coef_snps_by_gene, paste0("/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/", cell_type, "/coef_snps_by_gene_", start_pos, "to", end_pos, ".rds"))
saveRDS(pseudotime_by_gene, paste0("/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/", cell_type, "/pseudotime_by_gene_", start_pos, "to", end_pos, ".rds"))
saveRDS(interaction_by_gene, paste0("/work-zfs/abattle4/prashanthi/sc-endo/results/eQTL_calling/", cell_type, "/interaction_by_gene_", start_pos, "to", end_pos, ".rds"))

