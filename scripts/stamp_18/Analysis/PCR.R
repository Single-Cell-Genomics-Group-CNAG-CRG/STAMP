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

.pcr <- \(sce, x) {
  y <- reducedDim(sce, "PCA")
  z <- summary(lm(y ~ sce[[x]]))
  r2 <- sapply(z, \(.) .$adj.r.squared)
  data.frame(x, pc=seq_along(r2), r2)
}
# for multiple variables (mock code):
xs <- c("cell_area","sum","detected","sample")
df <- do.call(rbind, lapply(xs, \(x) .pcr(sce, x)))


pal_pcr <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_pcr) <- unique(df$x)

pcr <- ggplot(df, aes(x = pc, y = r2, color = x)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() + 
  scale_x_continuous(breaks = unique(df$pc)) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 25),
        aspect.ratio = 1/1.5) +
  scale_color_manual(values = pal_pcr) +
  scale_fill_manual(values = pal_pcr) +
  labs(x = "PC", color = "Metadata")

outdir <- glue("{plt_dir}/{stamp}/QC")
pdf(glue("{outdir}/pcr.pdf"), width = 8)
pcr
dev.off()
