rm(list=ls())
library(ggplot2)
plots.dir <- "../figures/"
egenes <- c()
for (cell_type in c("day1", "day3", "iPSC", "mesendo", "defendo")) {
    sigs <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/30pc/significant_hits.rds"))
    egenes <- c(egenes, nrow(sigs))
}
pdf(paste0(plots.dir, cell_type, "_egenes_v_type.pdf"))
data <- data.frame("types"=c("day1", "day3", "iPSC", "mesendo", "defendo"), "egenes"=egenes)
ggplot(data, aes(x=types, y=egenes)) + geom_col()
dev.off()
