rm(list = ls())

# This script merges cells which belong to the same individual, experiment and day
# Further we format the genotype matrix
# We also compute expression PCs to be used in downstream analysis for our samples

library(tidyverse)
library(Seurat)
library(genoscapeRtools)

inputArgs <-  commandArgs(TRUE)
cell_type <- inputArgs[1]

dat.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/data/"

expr <- read.csv(paste0(dat.dir, "expr/", cell_type, ".csv"))
rownames(expr) <- expr$X
expr$X <- NULL

metadata <- read.csv(paste0(dat.dir, "metadata/iPSC.csv"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$X
metadata$X <- NULL

geno.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/"

chr.matrices <- list()
for(i in c(1:22)){
  chr.matrices[[i]] <- read_012(paste0(geno.dir, "chr", as.character(i), "/combined.chr", as.character(i), ".common"), gz = FALSE)
}

geno_matrix <- t(chr.matrices[[1]])
for(i in c(2:22)){
  geno_matrix <- rbind(geno_matrix, t(chr.matrices[[i]]))
}

# Subset only to those individuals for whom we have genotype data
metadata <- metadata[metadata$donor_long_id %in% colnames(geno_matrix), ]
expr <- expr[rownames(expr) %in% rownames(metadata), ]
metadata <- metadata[match(rownames(metadata), rownames(expr)), ]
expr <- t(expr)

# Aggregate cells by donor, day of collection, experiment
metadata$agg_group <- paste(metadata$donor_long_id, metadata$experiment, metadata$day, sep = "_")
barplot(table(metadata$agg_group), names= FALSE, xlab = "Samples", ylab = "Number of cells", col = "light  gray", main = "iPSC cells")

colnames(expr) <- metadata$agg_group
uniq_samples <- unique(metadata$agg_group)
sample_count_df = sapply(uniq_samples, function(sample_id){
  sample_count_df = as.matrix(expr[ ,colnames(expr) == sample_id])
  sample_count_average_df = rowMeans(sample_count_df)
  return(sample_count_average_df)
})


saveRDS(sample_count_df, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr.rds"))

individuals <- c()
i <- 1
for(sample_id in uniq_samples){
  individuals[i] <- metadata$donor_long_id[metadata$agg_group == sample_id]
  i <- i + 1
}

geno_matrix_samples <- matrix(geno_matrix[ ,individuals[1]])
for(i in c(2:length(individuals))){
  geno_matrix_samples <- cbind(geno_matrix_samples, geno_matrix[ ,individuals[i]])
}

saveRDS(sample_count_df, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype.rds"))

# Compute expression PCs by approximating
expr.usv <- rsvd(scale(t(expr)), k = 100)
eigen_vectors <- expr.usv$u
rownames(eigen_vectors) <- colnames(expr)
eigen_vectors <- eigen_vectors[ ,1:10]
colnames(eigen_vectors) <- c("V1","V2","V3","V4","V5","V6","V7", "V8", "V9", "V10")

saveRDS(eigen_vectors, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr_PCs.rds"))

# Read the genotype PCs 
geno.pcs <- readRDS("/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/genotype_PCs.rds")
geno.pcs.augmented <- geno.pcs[individuals[1], ]
for(i in c(2:length(individuals))){
  geno.pcs.augmented <- rbind(geno.pcs.augmented, geno.pcs[individuals[i], ])
}
rownames(geno.pcs.augmented) <- uniq_samples

saveRDS(geno.pcs.augmented, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype_PCs.rds"))

# Map cells to individuals
random_effects_matrix <- matrix(0, nrow = length(uniq_samples), ncol = length(unique(individuals)))
rownames(random_effects_matrix) <- uniq_samples
colnames(random_effects_matrix) <- unique(individuals)
for(i in c(1:length(uniq_samples))){
  random_effects_matrix[uniq_samples[i], individuals[i]] <- random_effects_matrix[uniq_samples[i], individuals[i]]  + 1
}

saveRDS(random_effects_matrix, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/sample_assignments.rds"))


