rm(list=ls())

for cell_type in c("day1", "day3", "iPSC", "mesendo", "defendo") {
    for npcs in c(10, 15, 20, 30) {
	all <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/all_hits.rds"))
	sig <- readRDS(paste0("../results/eQTL_calling/", cell_type, "/", npcs, "pc/significant_hits.rds"))
	:
