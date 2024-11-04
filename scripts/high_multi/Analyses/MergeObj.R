# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(pheatmap)
library(SingleCellExperiment)
library(glue)
library(here)
library(qs)
library(ggvenn)

stamps <- c("stamp_17","stamp_18","stamp_15","stamp_16")
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# load data
s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/clust_sce.qs"), nthreads  = 8)
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"), nthreads = 8)
s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/PreProcNew.qs"), nthreads = 8)
s16 <- qread(glue("{proj_dir}/data/{stamps[4]}/processed/PreProcNew.qs"), nthreads = 8)


ft <- intersect(rownames(s15),rownames(s18))
s15 <- s15[ft,]
s16 <- s16[ft,]
s17 <- s17[ft,]
s18 <- s18[ft,]

# add replicate data
s15$replicate <- "s15"
s16$replicate <- "s16"
s17$replicate <- "s17"
s18$replicate <- "s18"

# add tech data
s15$tech <- "CosMx"
s16$tech <- "CosMx"
s17$tech <- "Xenium"
s18$tech <- "Xenium"

colnames(s15) <- paste0(colnames(s15),"_",s15$replicate,"_",s15$tech)
colnames(s16) <- paste0(colnames(s16),"_",s16$replicate,"_",s16$tech)
colnames(s17) <- paste0(colnames(s17),"_",s17$replicate,"_",s17$tech)
colnames(s18) <- paste0(colnames(s18),"_",s18$replicate,"_",s18$tech)

# Bind cosmx and xenium 
xenium <- cbind(s17,s18)
cosmx <- cbind(s15,s16)
# free memory
rm(s15)
rm(s16)
rm(s17)
rm(s18)
gc()

colData(xenium) <- colData(xenium)[,c("sample","sum","detected","replicate","tech")]
colData(cosmx) <- colData(cosmx)[,c("sample","sum","detected","replicate","tech")]

# create new clean obj
xenium <- SingleCellExperiment(assays = list(counts = counts(xenium)),
                               colData = colData(xenium))
cosmx <-  SingleCellExperiment(assays = list(counts = counts(cosmx)),
                               colData = colData(cosmx))

sce <- cbind(xenium,cosmx)

sce <- sce[,colnames(sce) %!in% colnames(sce)[duplicated(colnames(sce))]] # remove 30 duplicated cells
# free mem
rm(xenium)
rm(cosmx)
gc()

# Find sample names shared across CosMx and Xenium
shared <- names(which(tapply(sce$tech, sce$sample, function(x) length(unique(x)) > 1)))
# subset for those
sub <- sce[,sce$sample %in% shared]

# save
res_dir <- glue("{proj_dir}/data/high_multi/processed")
dir.create(res_dir,showWarnings = F,recursive = T)
qsave(sub, file = glue("{res_dir}/merged.qs"), nthreads = 8)
