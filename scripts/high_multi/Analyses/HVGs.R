# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

stamps <- c("stamp_17","stamp_18")
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/clust_sce.qs"))
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"))

s17 <- s17[,sample(colnames(s17),10000)]

hvg_s17 <- list()
for(s in unique(s17$sample)){
  
  sub <- s17[,s17$sample == s]
  logNormCounts(sub)
dec.var_s17 <- modelGeneVar(sub, BPPARAM = bp)
hvg_s17[[s]] <- getTopHVGs(dec.var_s17, fdr.threshold = 1)
hvg_s17[[s]] <- hvg_s17[[s]][1:2000]
}


s18 <- s18[,sample(colnames(s18),10000)]

hvg_s18 <- list()
for(s in unique(s18$sample)){
  
  sub <- s18[,s18$sample == s]
  logNormCounts(sub)
  dec.var_s18 <- modelGeneVar(sub, BPPARAM = bp)
  hvg_s18[[s]] <- getTopHVGs(dec.var_s18, fdr.threshold = 1)
  hvg_s18[[s]] <- hvg_s18[[s]][1:2000]
}

length(intersect(hvg_s17[["A"]],hvg_s18[["A"]])) / length(c(hvg_s17[["A"]],hvg_s18[["A"]]))


dec.var_s17 <- modelGeneVar(s17, block = s17$sample, BPPARAM = bp)
dec.var_s18 <- modelGeneVar(s18, block = s18$sample, BPPARAM = bp)

hvg_s17 <- getTopHVGs(dec.var_s17, fdr.threshold = 1)
hvg_s18 <- getTopHVGs(dec.var_s18, fdr.threshold = 1)

length(intersect(hvg_s17,hvg_s18)) / n_distinct(c(hvg_s17,hvg_s18))









