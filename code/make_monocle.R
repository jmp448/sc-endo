library(monocle3)
library(tibble)
library(readr)
library(Matrix)

expression_matrix <- read_tsv("../data/counts.corrected.tsv", col_names=T)
expression_matrix <- column_to_rownames(expression_matrix, var="X1")
expression_matrix <- as.matrix(expression_matrix)

cell_annotation <- read_tsv("../data/cell_metadata.corrected.tsv", col_names=T)
cell_annotation <- column_to_rownames(cell_annotation, var="X1")

gene_annotation <- data.frame(apply(data.frame(matrix(rownames(expression_matrix))), 1, substring, first=17), stringsAsFactors=F)
rownames(gene_annotation) <- rownames(expression_matrix)
colnames(gene_annotation) <- c("gene_short_name")

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_annotation,
                         gene_metadata = gene_annotation)

saveRDS(cds, "../data/monocle.rds")