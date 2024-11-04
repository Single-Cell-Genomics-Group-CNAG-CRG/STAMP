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
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"), nthreads = 8)
s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/PreProcNew.qs"), nthreads = 8)

ft <- intersect(rownames(s15),rownames(s18))

x <- list(
  CosMx = rownames(s15), 
  Xenium = rownames(s18)
)

venn <- ggvenn(
  x,
  fill_color = c("#b0cce2", "#f9b2aa"),
  stroke_size = 0.5,
  set_name_size = 4)

stamp <- "high_multi"
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/venn.pdf"),width = 4,height = 3)
venn
dev.off()

