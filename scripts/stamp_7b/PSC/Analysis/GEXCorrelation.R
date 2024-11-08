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
cosmx <- qread(glue("{res_dir}/PreProcNew.qs"))
cosmx$tech <- "cosmx"
# Flex data
sample <- "combined"
outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC/{sample}")
flex <- qread(glue("{outdir}/proc_sce.qs"), nthreads = 8)
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
colnames(fl) <- df[colnames(fl), "sample"]
#fl <- t(rowsum(t(fl), group = colnames(fl)))
fl <- as(fl, "dgCMatrix")
# Aggregate the cosmx matrix
df <- as.data.frame(colData(c))
colnames(cs) <- df[colnames(cs), "sample"]
#cs <- t(rowsum(t(cs), group = colnames(cs)))
cs <- as(cs, "dgCMatrix")

# function to aggregate and normalize
aggr_norm <- function(mat){
  res <- do.call(cbind, lapply(unique(colnames(mat)), function(col_name) {
    sub <- mat[, colnames(mat) %in% col_name, drop = FALSE]
    agg <- rowSums(sub) / ncol(sub)
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
  
  genes_iPSCs <- c("FGF2","POU5F1","WNT3","KDR")
  genes_meso <- c("GATA3","KRT19","FOXF1","BMP4","IGFBP7","TTN","EOMES")
  genes_ecto <- c("SNAI1","DLL1","SOX2","NRG1")

  # Select the appropriate marker genes based on the subset
  if (subset == "iPSC_parental") { feat_color <- genes_iPSCs }
  if (subset == "mesoderm") { feat_color <- genes_meso }
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
  
 if(subset == "iPSC_parental") {plt <- plt + scale_y_log10(limits = c(0.01,100), breaks = c(0.01,1,100))}
if(subset == "ectoderm" | subset == "mesoderm") {plt <- plt + scale_y_log10(limits = c(0.1,10), breaks = c(0.1,1,10))}
  
  return(plt)
}

# Combine the plots as before
gg_corr <- wrap_plots(
  corr_plot(df, "iPSC_parental"),
  corr_plot(df, "mesoderm"),
  corr_plot(df, "ectoderm"),
  nrow = 1
) +
  plot_layout(axis_titles = "collect")
#gg_corr


pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC_V2/gex_corr.pdf", width = 10, height = 3)
gg_corr
dev.off()
  

