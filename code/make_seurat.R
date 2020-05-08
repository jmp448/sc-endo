library(Seurat)
library(tibble)
library(readr)
library(Matrix)

expression_matrix <- read_tsv("../data/counts.corrected.tsv", col_names=T)
expression_matrix <- column_to_rownames(expression_matrix, var="X1")
rownames(expression_matrix) <- apply(matrix(rownames(expression_matrix)), 1, substring, first=17)

cell_annotation <- read_tsv("../data/cell_metadata.corrected.tsv", col_names=T)
cell_annotation <- column_to_rownames(cell_annotation, var="X1")

sobj <- CreateSeuratObject(expression_matrix,
                           meta.data = cell_annotation)

saveRDS(sobj, "../data/seurat.rds")