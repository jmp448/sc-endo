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
      OR <- log10((sig.in/(sig.tot-sig.in))/(bg.in/(bg.tot-bg.in)))
      p.fish <- fisher.test(matrix(c(sig.in, sig.tot-sig.in, bg.in, bg.tot-bg.in), nrow=2, ncol=2))$p.value
      if (p.fish <= 0.01) {
	lab="***"
      } else if (p.fish <= 0.05) {
        lab="**"
      } else if (p.fish <= 0.1) {
	lab="*"
      } else {
	lab=""
      }
      if (epi.mark=="H3K27me3" & context=="E018") {
        p.mat <- data.frame("marker"=epi.mark, "context"=context, "OR"=OR, "signif"=lab)
      } else {
        p.mat2 <- data.frame("marker"=epi.mark, "context"=context, "OR"=OR, "signif"=lab)
	p.mat <- rbind(p.mat, p.mat2)
      }
    }
  }
  png(paste0("../figures/", pseudotime, ".png"))
  p.mat$marker <- factor(p.mat$marker, levels=c("H3K4me1", "H3K4me3", "H3K27me3", "H3K36me3"))
  p.mat$context <- factor(p.mat$context, levels=contexts)
  plt <- ggplot(p.mat, aes(x=marker,
                    fill=context,
                    y=OR,
		    label=signif)) + 
    geom_bar(stat="identity", position="dodge") +
    xlab("Epigenetic Mark") + 
    ylab("Log10 Odds Ratio") +
    ggtitle(pseudotime) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(values=c("turquoise4", "burlywood2", "deeppink4")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
            panel.background = element_blank())
  plot(plt)
  dev.off()
}
