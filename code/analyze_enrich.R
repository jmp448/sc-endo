rm(list=ls())

library(R.utils)
library(ggplot2)

cell.types <- c("iPSC", "mesendo", "defendo")
ref.types <- c("E018", "E004", "E011")

for (epi.mark in c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3")) {
  for (i in c(1,2,3)) {
    sig.tot <- countLines(paste0("../results/eQTL_calling/", cell.types[i], "/10pc/significant_hits.bed"))
    bg.tot <- countLines(paste0("../results/eQTL_calling/", cell.types[i], "/10pc/background_hits.bed"))
    for (j in c(1,2,3)) {
      sig.in <- countLines(paste0("../results/func_analysis/", cell.types[i], "/", ref.types[j], "-", epi.mark, ".sig.bed"))
      bg.in <- countLines(paste0("../results/func_analysis/", cell.types[i], "/", ref.types[j], "-", epi.mark, ".bg.bed"))
      OR <- log10((sig.in/(sig.tot-sig.in))/(bg.in/(bg.tot-bg.in))) 
      fishy <- matrix(c(sig.in, sig.tot-sig.in, bg.in, bg.tot-bg.in), nrow=2, ncol=2)
      result <- fisher.test(fishy)
      print(paste0(epi.mark, " context ", ref.types[j], " cell type ", cell.types[i], " pval: ", result$p.value)) 
      if (i==1 & j==1) {
        p.mat <- data.frame("type"=cell.types[i], "context"=ref.types[j], "log10_OR"=OR)
      } else {
        p.mat2 <- data.frame("type"=cell.types[i], "context"=ref.types[j], "log10_OR"=OR)
        p.mat <- rbind(p.mat, p.mat2)
      }
    }
  }
  png(paste0("../figures/static.", epi.mark, ".png"))
  p.mat$type <- factor(p.mat$type, levels=c("iPSC", "mesendo", "defendo"))
  p.mat$context <- factor(p.mat$context, levels=c("E018", "E004", "E011"))
  plt <- ggplot(p.mat, aes(x=type,
                    y=context,
                    fill=log10_OR)) + 
    geom_raster()
  plot(plt)
  dev.off()
}
