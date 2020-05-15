rm(list = ls())

# This script merges cells which belong to the same individual, experiment and day
# Further we format the genotype matrix
# We also compute expression PCs to be used in downstream analysis for our samples

library(rsvd)
library(tidyverse)
library(Seurat)
library(genoscapeRtools)
require(data.table)

#inputArgs <-  commandArgs(TRUE)
#cell_type <- inputArgs[1]

cell_type <- "mesendo"

dat.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/data/"

expr <- fread(paste0(dat.dir, "expr/", cell_type, ".csv"), header = T)
#read.csv(paste0(dat.dir, "expr/", cell_type, ".csv"))
#rownames(expr) <- expr$V1
#expr$V1 <- NULL

metadata <- read.csv(paste0(dat.dir, "metadata/", cell_type, ".csv"), stringsAsFactors = FALSE)
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
expr <- expr[expr$V1 %in% rownames(metadata), ]
metadata <- metadata[match(rownames(metadata), expr$V1), ]
rownames(expr) <- expr$V1
expr$V1 <- NULL
expr <- t(expr)
colnames(expr) <- rownames(metadata)

# Aggregate cells by donor, day of collection, experiment
metadata$agg_group <- paste(metadata$donor_long_id, metadata$experiment, metadata$day, sep = "_")
barplot(table(metadata$agg_group), names= FALSE, xlab = "Samples", 
        ylab = "Number of cells", col = "light  gray", main = cell_type)

colnames(expr) <- metadata$agg_group
uniq_samples <- unique(metadata$agg_group)
sample_count_df = sapply(uniq_samples, function(sample_id){
  sample_count_df = as.matrix(expr[ ,colnames(expr) == sample_id])
  sample_count_average_df = rowMeans(sample_count_df)
  return(sample_count_average_df)
})




individuals <- c()
day <- c()
exp <- c()
i <- 1
for(sample_id in uniq_samples){
  individuals[i] <- metadata$donor_long_id[metadata$agg_group == sample_id]
  day[i] <- metadata$day[metadata$agg_group == sample_id]
  exp[i] <- metadata$experiment[metadata$agg_group == sample_id]
  i <- i + 1
}

geno_matrix_samples <- matrix(geno_matrix[ ,individuals[1]])
for(i in c(2:length(individuals))){
  geno_matrix_samples <- cbind(geno_matrix_samples, geno_matrix[ ,individuals[i]])
}
colnames(geno_matrix_samples) <- uniq_samples
saveRDS(geno_matrix_samples, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype.rds"))

# Compute expression PCs by approximating
scaled.expr <- scale(t(sample_count_df))
genes.na <- colnames(scaled.expr)[colSums(is.na(scaled.expr)) > 0]
scaled.expr <- scaled.expr[ ,!colnames(scaled.expr) %in% genes.na]
expr.usv <- svd(scaled.expr)
eigen_vectors <- expr.usv$u
pc.percent <- ((expr.usv$d^2)/sum((expr.usv$d^2)))*100

sample_count_df <- sample_count_df[!rownames(sample_count_df) %in% genes.na, ]
saveRDS(sample_count_df, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/expr.rds"))

plots.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/figures/"
pdf(paste0(plots.dir, paste0(cell_type, " Scree plot.pdf")))
plot(pc.percent, pch = 20, xlim = c(1,20), xlab = "PC Number", 
     ylab = "% Variance explained", main = paste0("Scree plot: ", cell_type))
lines(pc.percent)
dev.off()

if(cell_type == "all_cells"){
custom_col <- c("#009E73", "#0072B2", "#D55E00", "#CC79A7")
plot(expr.usv$u[ ,1], expr.usv$u[ ,2], col = custom_col[as.factor(day)], xlab = "PC1", ylab = "PC2", main = "All cells", pch = 20)
legend("bottomright", fill = custom_col, legend = c("Day0", "Day1", "Day2", "Day3"), cex = 0.7)
}

#plot(expr.usv$u[ ,1], expr.usv$u[ ,2], col = as.factor(day), xlab = "PC1", ylab = "PC2", main = "Defendo (day)", pch = 20)

rownames(eigen_vectors) <- colnames(sample_count_df)
eigen_vectors <- eigen_vectors[ ,1:35]
colnames(eigen_vectors) <- c("V1","V2","V3","V4","V5","V6","V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", 
                             "V16","V17","V18","V19","V20","V21","V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", 
                             "V31", "V32", "V33", "V34", "V35")
saveRDS(eigen_vectors, paste0("../data/eQTL_calling/", cell_type, "/expr_35PCs.rds"))

plots.dir <- "../figures/"
pdf(paste0(plots.dir, cell_type, "_35pcs.pdf"))
plot(pc.percent, pch = 20, xlim = c(1,35), xlab = "PC Number", 
     ylab = "% Variance explained", main = paste0("Scree plot: ", cell_type))
lines(pc.percent)
dev.off()


# Read the genotype PCs 
geno.pcs <- readRDS("/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/genotype_PCs.rds")
geno.pcs.augmented <- geno.pcs[individuals[1], ]
for(i in c(2:length(individuals))){
  geno.pcs.augmented <- rbind(geno.pcs.augmented, geno.pcs[individuals[i], ])
}
rownames(geno.pcs.augmented) <- uniq_samples

saveRDS(geno.pcs.augmented, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/genotype_PCs.rds"))

# Map cells to individuals
random_effects <- data.frame(uniq_samples, individuals)
colnames(random_effects) <- c("sample_id", "individuals")
saveRDS(random_effects, paste0("/work-zfs/abattle4/prashanthi/sc-endo/data/eQTL_calling/", cell_type, "/sample_assignments.rds"))


