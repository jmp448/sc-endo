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
      fishy <- matrix(c(sig.in, sig.tot-sig.in, bg.in, bg.tot-bg.in), nrow=2, ncol=2)
      pval <- fisher.test(fishy)$p.value
      odds.r <- (sig.in/(sig.tot-sig.in))/(bg.in/(bg.tot-bg.in))
      if (i==1 & j==1) {
        p.mat <- data.frame("type"=cell.types[i], "context"=ref.types[j], "enrichment"=odds.r)
      } else {
        p.mat2 <- data.frame("type"=cell.types[i], "context"=ref.types[j], "enrichment"=odds.r)
        p.mat <- rbind(p.mat, p.mat2)
      }
    }
  }
  png(paste0("../figures/static.", epi.mark, ".png"))
  plt <- ggplot(p.mat, aes(x=factor(type, levels=c("iPSC", "mesendo", "defendo")),
                    y=factor(context, levels=c("E018", "E004", "E011")),
                    fill=enrichment)) + 
    geom_raster()
  plot(plt)
  dev.off()
}
