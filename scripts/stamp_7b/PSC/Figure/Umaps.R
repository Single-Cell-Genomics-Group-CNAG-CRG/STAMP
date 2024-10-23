library(dplyr)
library(patchwork)
library(ggplot2)
library(qs)
library(glue)
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
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~tech, ncol = 1)

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/comb_umaps.pdf", width = 4, height = 4)
umap_clusters
dev.off()


# Proportions and cell numbers
prop_plt <- \(samp){
    
  cd_f <- as.data.frame(colData(flex[,flex$sample == samp]))
  cd_c <- as.data.frame(colData(cosmx[,cosmx$sample == samp]))
  
  cd_c <- cd_c %>% select(tech,cluster, sample)
  cd_f <- cd_f %>% select(tech,cluster, sample)
  
  cd <- rbind(cd_c,cd_f)
  
  df <- as.data.frame(table(cd$tech,cd$cluster))
  df <- df %>% group_by(Var1) %>% mutate(perc = round(Freq/sum(Freq),2)) %>% ungroup()
  gg_prop <- ggplot(df, aes(x = Var1, y = perc, fill = Var2)) + 
    geom_col() + 
    scale_fill_manual(values = pal) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1,
                                     vjust = 1, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          title = element_text(size = 18, color = "black"),
          aspect.ratio = 3/1,
          legend.position = "none") +
    labs(x = "Tech", y = "Proportions", fill = "Cluster", title = samp) 
  return(gg_prop)
}

gg_prop <- wrap_plots(prop_plt("mesoderm"), 
           prop_plt("ectoderm"), 
           prop_plt("iPSC_parental"), 
           nrow = 2,ncol = 2) +
  plot_layout(guides = "collect", axis_titles = "collect")



pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/prop.pdf", width = 8, height = 8)
gg_prop
dev.off()

## qcmetrics
cd <- as.data.frame(colData(cosmx))
cd$cluster <- factor(cd$cluster, levels = names(sort(tapply(cd$sum, cd$cluster, mean), decreasing = FALSE)))

# Plot
qcmet <- \(var){
  if(var == "sum") {title = "nCount"}
  if(var == "detected") {title = "nFeature"}
  if(var == "Area.um2") {title = "Area.um2"}
  
  plt <- ggplot(cd, aes(x = cluster, y = !!sym(var),fill = cluster,)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.1) +
    scale_y_log10() + 
    scale_fill_manual(values = pal) + 
    scale_color_manual(values = pal) +
    labs(y = title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 15),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          aspect.ratio = 1/2)
  
  if(var == "sum") {plt <- plt + scale_y_log10()}
  if(var == "detected") {plt <- plt + scale_y_log10()}
  return(plt)
}

wrap_plots(
  qcmet("sum"),
  qcmet("detected"),
  qcmet("Area.um2"),
  ncol =1) +
  plot_layout(guides = "collect", axis_titles = "collect")
