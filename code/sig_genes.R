rm(list=ls())

plots.dir <- "../figures/"
egenes <- c()
for (cell_type in c("day1", "day3", "iPSC", "mesendo", "defendo")) {
    sigs <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/30pc/significant_hits.rds"))
    egenes <- c(egenes, nrow(sigs))
}
pdf(paste0(plots.dir, cell_type, "_egenes_v_type.pdf"))
data <- table(egenes, rownames=c("day1", "day3", "iPSC", "mesendo", "defendo"))
barplot(data, pch = 20, xlab = "Cell Type", ylab = "Significant eGenes")
dev.off()
