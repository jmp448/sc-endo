rm(list=ls())

library(R.utils)
library(ggplot2)
library(stats)

contexts <- c("E018", "E004", "E011")

for (pseudotime in c("stegle", "UMAP")) {
  sig.tot <- countLines(paste0("../results/eQTL_calling/all_cells/", pseudotime, "/10pc/significant_hits.bed"))
  bg.tot <- countLines(paste0("../results/func_analysis/static.merged.bed"))
  for (epi.mark in c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3")) {
    for (context in contexts) {
      sig.in <- countLines(paste0("../results/func_analysis/all_cells/", pseudotime, "/", context, "-", epi.mark, ".sig.bed"))
      bg.in <- countLines(paste0("../results/func_analysis/", context, "-", epi.mark, ".static.bed"))
      OR <- (sig.in/(sig.tot-sig.in))/(bg.in/(bg.tot-bg.in))
      if (epi.mark=="H3K27me3" & context=="E018") {
        p.mat <- data.frame("marker"=epi.mark, "context"=context, "OR"=OR)
      } else {
        p.mat2 <- data.frame("marker"=epi.mark, "context"=context, "OR"=OR)
	p.mat <- rbind(p.mat, p.mat2)
      }
    }
  }
  png(paste0("../figures/", pseudotime, ".png"))
  p.mat$marker <- factor(p.mat$marker, levels=c("H3K4me1", "H3K4me3", "H3K27me3", "H3K36me3"))
  p.mat$context <- factor(p.mat$context, levels=contexts)
  plt <- ggplot(p.mat, aes(x=marker,
                    fill=context,
                    y=OR)) + 
    geom_bar(stat="identity", position="dodge") +
    xlab("Epigenetic Mark") + 
    ylab("Odds Ratio") +
    ggtitle(pseudotime) +
    geom_hline(yintercept=1, linetype="dashed")
  plot(plt)
  dev.off()
}
