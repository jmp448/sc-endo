rm(list = ls())

# This script merges cells which belong to the same individual, experiment and day
# Further we format the genotype matrix
# We also compute expression PCs to be used in downstream analysis for our samples

library(rsvd)
library(tidyverse)
library(Seurat)
require(data.table)

#inputArgs <-  commandArgs(TRUE)
#cell_type <- inputArgs[1]

cell_types <- c("day1", "day3")

for (cell_type in cell_types) {
  sample_count_df <- readRDS(paste0("../data/eQTL_calling/", cell_type, "/expr.rds"))
  # Compute expression PCs by approximating
  expr.usv <- rsvd(scale(t(sample_count_df)))
  eigen_vectors <- expr.usv$u
  pc.percent <- ((expr.usv$d^2)/sum((expr.usv$d^2)))*100
  
  plots.dir <- "../figures/"
  pdf(paste0(plots.dir, cell_type, "_35pcs.pdf"))
  plot(pc.percent, pch = 20, xlim = c(1,35), xlab = "PC Number", 
       ylab = "% Variance explained", main = paste0("Scree plot: ", cell_type))
  lines(pc.percent)
  dev.off()
  
  rownames(eigen_vectors) <- colnames(sample_count_df)
  eigen_vectors <- eigen_vectors[ ,1:35]
  colnames(eigen_vectors) <- c("V1","V2","V3","V4","V5","V6","V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35")
  saveRDS(eigen_vectors, paste0("../data/eQTL_calling/", cell_type, "/expr_35PCs.rds"))
}
