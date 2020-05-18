rm(list = ls())

library(rsvd)
library(tidyverse)
library(Seurat)
require(data.table)
library(RColorBrewer)
library(genoscapeRtools)
library(nlme)
library(pbkrtest)

dat.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/data/"

metadata <- read.csv(paste0(dat.dir, "metadata.pseudotime.tsv"), stringsAsFactors = FALSE, sep = "\t")
metadata$agg_id <- paste(metadata$donor_short_id, metadata$experiment, sep = "-")

uniq.ids <- unique(metadata$agg_id)
pseudotime_0 <- c()
pseudotime_1 <- c()
pseudotime_2 <- c()
pseudotime_3 <- c()

for(id in uniq.ids){
  pseudotime_0 <- c(pseudotime_0, mean(metadata$pseudo[metadata$agg_id == id & metadata$day == "day0"]))
  pseudotime_1 <- c(pseudotime_1, mean(metadata$pseudo[metadata$agg_id == id & metadata$day == "day1"]))
  pseudotime_2 <- c(pseudotime_2, mean(metadata$pseudo[metadata$agg_id == id & metadata$day == "day2"]))
  pseudotime_3 <- c(pseudotime_3, mean(metadata$pseudo[metadata$agg_id == id & metadata$day == "day3"]))
}

df <- data.frame(uniq.ids, pseudotime_0, pseudotime_1, pseudotime_2, pseudotime_3)
count_na <- rowSums(is.na(df[ ,2:5]))
df <- df[count_na == 0, ]
rybPal <- colorRampPalette(brewer.pal(n = 11, name = 'RdYlBu'))

df$col <- rybPal(10)[as.numeric(cut(df$pseudotime_3,breaks = 10))]

pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/differentiation_efficiency.pdf")
plot(c(0,1,2,3), df[1, 2:5], col = df$col[1], pch = 20, xlab = "Time point", 
     ylab = "Pseudotime", type = 'l', xaxt = "n", ylim = c(0,0.8), main = "Differentiation efficiency")
axis(1, at=0:3, labels=c("Day0", "Day1", "Day2", "Day3"))
for(i in c(2:dim(df)[1])){
  lines(c(0,1,2,3), df[i, 2:5], col = df$col[i])
}
legend("bottomright",title="Decile",legend=c(1:10),col =rybPal(10),pch=15, ncol = 2)
dev.off()

# Get genotype
df$uniq.ids <- as.character(df$uniq.ids)
df$individual <- gsub("-.*", "", df$uniq.ids)
df$exp <- gsub(".*-", "", df$uniq.ids)

significant_hits_ipsc <- readRDS("~/data/prashanthi/sc-endo/results/eQTL_calling/iPSC/10pc/significant_hits.rds")
significant_hits_mesendo <- readRDS("~/data/prashanthi/sc-endo/results/eQTL_calling/mesendo/10pc/significant_hits.rds")
significant_hits_defendo <- readRDS("~/data/prashanthi/sc-endo/results/eQTL_calling/defendo/10pc/significant_hits.rds")

geno.dir <- "/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/"

chr.matrices <- list()
for(i in c(1:22)){
  chr.matrices[[i]] <- read_012(paste0(geno.dir, "chr", as.character(i), "/combined.chr", as.character(i), ".common"), gz = FALSE)
}

geno_matrix <- t(chr.matrices[[1]])
for(i in c(2:22)){
  geno_matrix <- rbind(geno_matrix, t(chr.matrices[[i]]))
}

geno_matrix <- t(geno_matrix)
geno_present <- c()
for(ind in df$individual){
  geno_present <- c(geno_present, length(grep(ind, rownames(geno_matrix))))
}
df <- df[!geno_present == 0, ]

snps_iPSC <- c()
pvalue_iPSC <- c()

significant_hits_ipsc$snps <- as.character(significant_hits_ipsc$snps)
significant_hits_ipsc <- significant_hits_ipsc[!significant_hits_ipsc$snps == "geno_subset", ]

for(isnp in significant_hits_ipsc$snps){
  isnp <- gsub("`", "", isnp)
  isnp <- gsub("`", "", isnp)
  igeno_matrix <- geno_matrix[ ,isnp]
  igeno <- c()
  for(ind in df$individual){
    igeno <- c(igeno, igeno_matrix[grep(ind, names(igeno_matrix))]) }
  df_lm <- cbind(df, igeno)
  geno.LMM <- lmer(pseudotime_3 ~ igeno + (1 |individual) + (1| exp) , data = df_lm)
  icoefs <- data.frame(coef(summary(geno.LMM)))
  
  df.KR <- get_Lb_ddf(geno.LMM, fixef(geno.LMM))
  icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
  snps_iPSC <- c(snps_iPSC, isnp)
  pvalue_iPSC <- c(pvalue_iPSC, icoefs$p.KR[2])
}

iPSC <- data.frame(snps_iPSC, pvalue_iPSC)
iPSC$adj_pvalue <-p.adjust(iPSC$pvalue_iPSC, method = "BH")


snps_mesendo <- c()
pvalue_mesendo <- c()

significant_hits_mesendo$snps <- as.character(significant_hits_mesendo$snps)
significant_hits_mesendo <- significant_hits_mesendo[!significant_hits_mesendo$snps == "geno_subset", ]

for(isnp in significant_hits_mesendo$snps){
  isnp <- gsub("`", "", isnp)
  isnp <- gsub("`", "", isnp)
  igeno_matrix <- geno_matrix[ ,isnp]
  igeno <- c()
  for(ind in df$individual){
    igeno <- c(igeno, igeno_matrix[grep(ind, names(igeno_matrix))]) }
  df_lm <- cbind(df, igeno)
  geno.LMM <- lmer(pseudotime_3 ~ igeno + (1 |individual) + (1| exp) , data = df_lm)
  icoefs <- data.frame(coef(summary(geno.LMM)))
  
  df.KR <- get_Lb_ddf(geno.LMM, fixef(geno.LMM))
  icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
  snps_mesendo <- c(snps_mesendo, isnp)
  pvalue_mesendo <- c(pvalue_mesendo, icoefs$p.KR[2])
}

mesendo <- data.frame(snps_mesendo, pvalue_mesendo)
mesendo$adj_pvalue <-p.adjust(mesendo$pvalue_mesendo, method = "BH")

snps_defendo <- c()
pvalue_defendo <- c()

significant_hits_defendo$snps <- as.character(significant_hits_defendo$snps)
significant_hits_defendo <- significant_hits_defendo[!significant_hits_defendo$snps == "geno_subset", ]

for(isnp in significant_hits_defendo$snps){
  isnp <- gsub("`", "", isnp)
  isnp <- gsub("`", "", isnp)
  igeno_matrix <- geno_matrix[ ,isnp]
  igeno <- c()
  for(ind in df$individual){
    igeno <- c(igeno, igeno_matrix[grep(ind, names(igeno_matrix))]) }
  df_lm <- cbind(df, igeno)
  geno.LMM <- lmer(pseudotime_3 ~ igeno + (1 |individual) + (1| exp) , data = df_lm)
  icoefs <- data.frame(coef(summary(geno.LMM)))
  
  df.KR <- get_Lb_ddf(geno.LMM, fixef(geno.LMM))
  icoefs$p.KR <- 2 * (1 - pt(abs(icoefs$t.value), df.KR))
  snps_defendo <- c(snps_defendo, isnp)
  pvalue_defendo <- c(pvalue_defendo, icoefs$p.KR[2])
}

defendo <- data.frame(snps_defendo, pvalue_defendo)
defendo$adj_pvalue <-p.adjust(defendo$pvalue_defendo, method = "BH")

iPSC$snps_iPSC <- gsub("--", ".", iPSC$snps_iPSC)
iPSC$snps_iPSC <- as.numeric(iPSC$snps_iPSC)
mesendo$snps_mesendo <- gsub("--", ".", mesendo$snps_mesendo)
mesendo$snps_mesendo <- as.numeric(mesendo$snps_mesendo)
defendo$snps_defendo <- gsub("--", ".", defendo$snps_defendo)
defendo$snps_defendo <- as.numeric(defendo$snps_defendo)


pdf("/work-zfs/abattle4/prashanthi/sc-endo/figures/Predicting_efficiency.pdf")
plot(iPSC$snps_iPSC , -log10(iPSC$adj_pvalue), col = "turquoise4",
pch = 20, axes = F, xlab = "Chromosome", ylab = "-log10(p value)", main = "Predicting diff. efficiency", ylim = c(0.05, 1.5))
points(mesendo$snps_mesendo, -log10(mesendo$adj_pvalue), col = "tan1", pch = 20)
points(defendo$snps_defendo, -log10(defendo$adj_pvalue), col = "deeppink4", pch = 20)
text(3.504122, -log10(min(defendo$adj_pvalue)), labels = "BHLHE40", pos = 4)
axis(1)
axis(2)
abline(h = 1.0, col = "black", lty = 2)
legend("topright", legend = c("iPSC", "mesendo", "defendo"), fill = c("turquoise4", "tan1", "deeppink4"))
dev.off()

