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


##  marker genes
mrk <- c("GATA3","IGFBP7","KRT80","KRT19",
         "DLL1","SOX2","BMP4",
         "WNT5A","PDGFRA","KDR","TTN","FOXF1","SNAI1",
         "WNT3","POU5F1","FGF2")

### CORR SPLIT BY CELL TYPE
# Take matrices
c <- cosmx
f <- flex
cs <- as(counts(c),"dgCMatrix")
fl <- as(counts(f), "dgCMatrix")

feat <- intersect(rownames(cs),rownames(fl)) # take features intersection
fl <- fl[feat, ] # subset 
cs <- cs[feat,] # subset

# Aggregate the flex  matrix
df <- as.data.frame(colData(f))
colnames(fl) <- df[colnames(fl), "cluster"]
#fl <- t(rowsum(t(fl), group = colnames(fl)))
fl <- as(fl, "dgCMatrix")
# Aggregate the cosmx matrix
df <- as.data.frame(colData(c))
colnames(cs) <- df[colnames(cs), "cluster"]
#cs <- t(rowsum(t(cs), group = colnames(cs)))
cs <- as(cs, "dgCMatrix")

# function to aggregate and normalize
aggr_norm <- function(mat){
  res <- do.call(cbind, lapply(unique(colnames(mat)), function(col_name) {
    sub <- mat[, colnames(mat) %in% col_name, drop = FALSE]
    agg <- rowSums(sub) / ncol(sub)
  #  agg <- rowSums(sub)
    agg
  }))
  colnames(res) <- unique(colnames(mat))
  return(res)
}

agg_fl <- as.data.frame(aggr_norm(fl))
agg_fl$gene <- rownames(agg_fl)
agg_cs <- as.data.frame(aggr_norm(cs))
agg_cs$gene <- rownames(agg_cs)


agg_cs <- melt(agg_cs, id.vars = "gene")
names(agg_cs)[names(agg_cs) == "value"] <- "cosmx"

agg_fl <- melt(agg_fl, id.vars = "gene")
names(agg_fl)[names(agg_fl) == "value"] <- "flex"


df <- merge(agg_cs, agg_fl, by = c("gene", "variable"), all = TRUE)


corr_plot <- function(df, subset) {
  
  genes_amnion_like <- c("GATA3", "IGFBP7", "KRT80", "KRT19")
  genes_BMP <- c("GATA3", "DLL1", "BMP4")
  genes_ecto <- c("NRG1", "SOX2")
  genes_late_meso <- c("BMP4", "WNT5A", "PDGFRA", "KDR", "TTN", "FOXF1", "SNAI1", "WNT3")
  genes_pluri <- c("WNT3", "POU5F1", "FGF2")

  # Select the appropriate marker genes based on the subset
  if (subset == "amnion-like") { feat_color <- genes_amnion_like }
  if (subset == "BMP-induced prog.") { feat_color <- genes_BMP }
  if (subset == "late meso.") { feat_color <- genes_late_meso }
  if (subset == "pluripotent") { feat_color <- genes_pluri }
  if (subset == "ectoderm") { feat_color <- genes_ecto }
  
  # Filter the data for the specific subset
  sub <- df[df$variable == subset, ]
  # Create a new variable to indicate if a gene is a marker
  sub$Is_Marker <- ifelse(sub$gene %in% feat_color, "Yes", "No")
  
  # Separate the data into marker and non-marker genes
  sub_no <- sub[sub$Is_Marker == "No", ]
  sub_yes <- sub[sub$Is_Marker == "Yes", ]
  
  # Plot
  plt <- ggplot(sub, aes(x = flex, y = cosmx)) +
    # Plot non-marker genes first
    ggrastr::rasterise(
      geom_point(
        data = sub_no, 
        aes(color = Is_Marker),
        shape = 16, 
        size = 2, 
        alpha = 0.8
      ), 
      dpi = 800
    ) +
    # Add correlation coefficient
    ggpubr::stat_cor(
      method = "spearman", 
      label.x.npc = "left", 
      label.y.npc = "top", 
      size = 3
    ) + 
    # Plot marker genes on top
    geom_point(
      data = sub_yes, 
      aes(color = Is_Marker), 
      shape = 16, 
      size = 2, 
      alpha = 0.8
    ) +
    # Add gene labels for marker genes with black text
    ggrepel::geom_text_repel(
      data = sub_yes, 
      aes(label = gene),
      size = 3, 
      color = "black",  # Changed from "red" to "black"
      max.overlaps = Inf
    ) +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(
      slope = 1, 
      intercept = 0, 
      color = "red4", 
      linetype = "dashed"
    ) +
    labs(
      title = subset, 
      color = "Is Marker",
      x = "Mean counts/gene - Flex", 
      y = "Mean counts/gene - CosMx"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 14, color = "black"), 
    #  plot.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"  # Remove legend for combining later
    ) +
    scale_color_manual(values = c("Yes" = "red", "No" = "grey"))
  
  return(plt)
}

# Combine the plots as before
gg_corr <- wrap_plots(
  corr_plot(df, "pluripotent"),
  corr_plot(df, "amnion-like"),
  corr_plot(df, "late meso."),
  corr_plot(df, "ectoderm"),
  nrow = 1
) +
  plot_layout(axis_titles = "collect")
#gg_corr


pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/gex_corr.pdf", width = 12, height = 3)
gg_corr
dev.off()


