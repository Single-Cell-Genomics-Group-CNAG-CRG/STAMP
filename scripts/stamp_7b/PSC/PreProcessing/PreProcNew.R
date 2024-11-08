suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(qs)
  library(parallel)
  library(scran)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
})

stamp <- "stamp_7b"
sample <- "iPSCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
sce <- qread(glue("{res_dir}/qc_sce.qs"))

sce

# Subset for panel intersection of CosMx and Flex
#panel_int <- qread(glue("{proj_dir}/data/{stamp}/processed/panel_int.qs"))
#sce <- sce[panel_int,]

# LogNorm
sce <- logNormCounts(sce)

# PCA
set.seed(101001)
sce <- fixedPCA(sce, subset.row = NULL)

num_pcs_to_retain <- 10
percent.var <- attr(reducedDim(sce), "percentVar")

# Create a data frame for ggplot
data <- data.frame(PC = 1:length(percent.var), Variance = percent.var)
# Plot
gg_var <- ggplot(data, aes(x = PC, y = Variance)) +
  geom_point() +
  xlab("PC") +
  ylab("Variance explained (%)") +
  geom_vline(xintercept = num_pcs_to_retain, color = "red") +
  theme_bw()
gg_var

reducedDim(sce, "PCA") <-  reducedDim(sce, "PCA")[,1:num_pcs_to_retain]
wh(6,5)
gg_pca <- plotPCA(sce, scattermore = TRUE, point_size = 2) + ggtitle("PCA")
gg_pca

## Run UMAP
set.seed(123)
sce <- runUMAP(sce, dimred="PCA", BPPARAM = bp)
gg_um <- plotReducedDim(sce, "UMAP", scattermore = TRUE, point_size = 2) 
gg_um

sce$sample <- factor(sce$sample, levels = c("iPSC_parental","endoderm","mesoderm","ectoderm"))
# Save plots
plots_dir <- glue("{plt_dir}/{stamp}/{sample}")
dir.create(plots_dir, showWarnings = F)
pdf(glue("{plots_dir}/PreProc.pdf"), width = 7, height = 5)
wrap_plots(gg_var,gg_pca,gg_um, ncol = 2)
dev.off()

# Save data
qsave(sce, file = glue("{res_dir}/PreProcNew.qs"))