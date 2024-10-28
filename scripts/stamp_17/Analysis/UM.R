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

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")


# Take metadata and umap coordinates
df <- as.data.frame(colData(sce))
um <- as.data.frame(reducedDim(sce,"UMAP"))
df <- cbind(df,um) %>% select(sample,label, UMAP1,UMAP2)
df <- df[sample(rownames(df)),]

# Plot
um <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample)) +
  # geom_point(shape = 16, size = 0.01) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1200) +
  scale_color_manual(values = pal) +
  labs(x = "UMAP1", y = "UMAP2", color = "Sample") +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text =element_blank(),
    #axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none")
 
outdir <- glue("{plt_dir}/{stamp}/QC")
pdf(glue("{outdir}/um.pdf"), height = 3, width = 3)
um
dev.off()
