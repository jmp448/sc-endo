rm(list=ls())

plots.dir <- "../figures/"
for (cell_type in c("day1", "day3", "iPSC", "mesendo", "defendo")) {
    pcs <- c(10,20,30)
    egenes <- c()
    for (npcs in pcs) {
	sigs <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/significant_hits.rds"))
	egenes <- c(egenes, nrow(sigs))
    }
    pdf(paste0(plots.dir, cell_type, "_egenes_v_pc.pdf"))
    plot(pcs, egenes, pch = 20, xlab = "Number of Expression PCs", ylab = "Significant eGenes")
    dev.off()
}
