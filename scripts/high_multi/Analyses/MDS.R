# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(pheatmap)
library(SingleCellExperiment)
library(glue)
library(here)
library(qs)


stamp <- "high_multi"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
res_dir <- glue("{proj_dir}/data/high_multi/processed")

sce <- qread(glue("{res_dir}/integrated_25pct_merged.qs"), nthreads = 8)

set.seed(123)
pal <- Polychrome::createPalette(31, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V",
                "W","X","Y","MX1","Q","R",
                "s15","s16","s17","s18")
pal[28] <- "#A6CEE3"
pal[29] <- "#1F78B4" 
pal[30] <- "#B2DF8A" 
pal[31] <- "#33A02C"


agg <- aggregateAcrossCells(sce, ids = sce$s_r_t, use.assay = "counts", statistics = "mean")

agg <- logNormCounts(agg)
agg <- runMDS(agg)

plotMDS(agg, color_by = "tech", point_size = 2) 
 # scale_color_manual(values = pal)

mds <- as.data.frame(reducedDim(agg,"MDS"))
mds$sample <- rownames(mds)

ggplot(mds, aes(x = V1, y = V2, color = technology)) +
  geom_point(size = 3) +
  scale_color_manual(values = pal) +
  coord_equal() +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(color = "black", size = 15),
        panel.grid = element_blank())



outdir <- glue("{plt_dir}/{stamp}/QC")
pdf(glue("{outdir}/mds.pdf"), height = 3, width = 3)
mds
dev.off()








