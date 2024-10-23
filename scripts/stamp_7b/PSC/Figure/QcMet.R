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

cd <- as.data.frame(colData(cosmx))
pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(df$cluster)
# Calculate mean of sum per cluster and arrange clusters in decreasing order
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
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
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

  
