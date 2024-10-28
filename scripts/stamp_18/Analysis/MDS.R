suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(scuttle)
  library(glue)
  library(qs)
})

# Data loading
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_17"
res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/clust_sce.qs"))
sce

sce <- sce[,sample(colnames(sce),1000)]
sce <- logNormCounts(sce, BPPARAM = bp)


agg <- aggregateAcrossCells(sce, ids = sce$sample, use.assay.type = "logcounts", statistics = "mean")

agg <- runMDS(agg)

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

mds <- as.data.frame(reducedDim(agg,"MDS"))
mds$sample <- rownames(mds)

mds <- ggplot(mds, aes(x = V1, y = V2, color = sample)) +
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








