library(dplyr)
library(patchwork)
library(ggplot2)
library(qs)
library(here)

stamp <- "stamp_7b"
sample <- "iPSCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# CosMx data
res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
cosmx <- qread(glue("{res_dir}/anno_sce.qs"))
cosmx$tech <- "cosmx"
# Flex data
sample <- "combined"
outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC/{sample}")
flex <- qread(glue("{outdir}/anno_sce.qs"), nthreads = 8)
flex$tech <- "flex"

# Take metadata and umap coordinates
df_c <- as.data.frame(colData(cosmx))
um_c <- as.data.frame(reducedDim(cosmx,"UMAP"))
df_c <- cbind(df_c,um_c) %>% select(sample,cluster, V1,V2,tech)
colnames(df_c)[colnames(df_c) == "V1"] <- "UMAP1"
colnames(df_c)[colnames(df_c) == "V2"] <- "UMAP2"
df_c <- df_c[sample(rownames(df_c)),]

df_f <- as.data.frame(colData(flex))
um_f <- as.data.frame(reducedDim(flex,"UMAP"))
df_f <- cbind(df_f,um_f) %>% select(sample,cluster, UMAP1,UMAP2,tech)
df_f <- df_f[sample(rownames(df_f)),]

df <- rbind(df_c,df_f)

# Plot

## Palette
pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(df$cluster)

umap_clusters <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 800) +
  scale_color_manual(values = pal) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster") +
  theme_bw() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.text =element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~tech)

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/comb_umaps.pdf", width = 6, height = 3)
umap_clusters
dev.off()
  
