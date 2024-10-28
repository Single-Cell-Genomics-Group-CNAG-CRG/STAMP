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

stamp <- "stamp_15"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/qc_sce.qs"))
sce

#sce <- sce[,sample(colnames(sce),10000)]
# logNormCounts
sce <- logNormCounts(sce)

# Feature selection
set.seed(0010101)
dec.var <- modelGeneVar(sce, BPPARAM = bp) # model gene var
hvg <- getTopHVGs(dec.var,fdr.threshold = 1) # select hvg on fdr
dec.var$hvg <- "no" # Assign to dec.var column for plot
dec.var$hvg[rownames(dec.var) %in% hvg] <- "yes"
gg_hvg <- plot_hvg(dec.var = dec.var, sub = stamp) # plot
gg_hvg

# PCA
set.seed(101001)
sce <- fixedPCA(sce, subset.row = hvg)

num_pcs_to_retain <- 20
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

# run UMAP
set.seed(123)
sce <- runUMAP(sce, dimred="PCA", BPPARAM = bp)
gg_um <- plotReducedDim(sce, "UMAP", scattermore = TRUE, point_size = 2) 
gg_um

# Save plots
plots_dir <- glue("{plt_dir}/{stamp}/processed")
dir.create(plots_dir, showWarnings = F,recursive = T)
pdf(glue("{plots_dir}/PreProc.pdf"), width = 7, height = 5)
wrap_plots(gg_var,gg_hvg, gg_pca, gg_um, ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A") + 
  plot_annotation(title = glue("{stamp}"), subtitle = glue("N = {scales::comma(ncol(sce))} cells"))
dev.off()

# Save data
qsave(sce, file = glue("{res_dir}/PreProcNew.qs"), nthreads = 8)
