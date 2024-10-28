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
  library(reshape2)
})

stamp <- "stamp_18"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/clust_sce.qs"))

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

df <- as.data.frame(colData(sce))
pc <- as.data.frame(reducedDim(sce,"PCA")) %>% select(PC1,PC2,PC3,PC4)
df <- cbind(df, pc)
df <- df[sample(rownames(df)),]

pca <- ggplot(df, aes(x = PC1, y = PC2, color = sample)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1200) +
  scale_color_manual(values = pal) +
  labs(x = "PC1", y = "PC2", color = "Sample") +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text =element_blank(),
    #axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none") +
  coord_equal()

outdir <- glue("{plt_dir}/{stamp}/QC")
pdf(glue("{outdir}/pca_plt.pdf"), height = 4, width = 4)
pca
dev.off()




