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
sce <- qread(glue("{res_dir}/clustered_sce.qs"))

um <- plotReducedDim(sce,"UMAP",color_by = "sample", scattermore = T) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/um_cosmx.pdf", width = 3, height = 2)
um
dev.off()


df <- as.data.frame(colData(sce))
um <- as.data.frame(reducedDim(sce,"UMAP"))
df <- cbind(df,um)
df <- df[sample(rownames(df)),]

# clusters umap
clusters_um <- ggplot(df, aes(x = V1, y = V2, color = label)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1000) +
 # scale_color_manual(values = pal) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster") +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text =element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# samples umap
samples_um <- ggplot(df, aes(x = V1, y = V2, color = sample)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1000) +
  # scale_color_manual(values = pal) +
  labs(x = "UMAP1", y = "UMAP2", color = "Sample") +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text =element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

#samples_um

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/um_cosmx.pdf", width = 6, height = 4)
samples_um
dev.off()

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/um_cosmx_clusters.pdf", width = 6, height = 4)
clusters_um
dev.off()
