# dependencies
library(ggplot2)
library(dplyr)
library(ggsignif)
library(pheatmap)
library(SingleCellExperiment)
library(glue)
library(here)
library(qs)
library(scater)
library(scuttle)
library(scran)
library(BiocSingular)

stamp <- "high_multi"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
res_dir <- glue("{proj_dir}/data/high_multi/processed")

sce <- qread(glue("{res_dir}/integrated_25pct_merged.qs"), nthreads = 8)

cline_info <- readxl::read_xlsx("/Users/emanuelepitino/Desktop/high_multi/CellLineandMarkersSTAMPs.xlsx")

cline_info <- cline_info[cline_info$`Letter ID` %in% c(unique(sce$sample)),]

# extract markers in named vector
c <- lapply(cline_info$`Potential Markers`, function(x) trimws(unlist(strsplit(x, ","))))

c <- unlist(c)
### CORR SPLIT BY CELL TYPE
# Take matrix
cs <- as(counts(sce),"dgCMatrix")

# Aggregate the matrix by tech
df <- as.data.frame(colData(sce))
colnames(cs) <- df[colnames(cs), "tech"]


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

agg_cs <- as.data.frame(aggr_norm(cs))
agg_cs$gene <- rownames(agg_cs)

agg_cs$Is_Marker <- ifelse(agg_cs$gene %in% c, "Yes", "No")


sub <- agg_cs

sub_no <- sub[sub$Is_Marker == "No", ]
sub_yes <- sub[sub$Is_Marker == "Yes", ]

# Plot
corr <- ggplot(sub, aes(x = Xenium, y = CosMx)) +
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
  geom_smooth(method = "lm", color ="red4",linetype = "dashed") +
  labs(
    title = "subset", 
    color = "Is Marker",
    x = "Mean counts/gene - Xenium", 
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
  scale_color_manual(values = c("Yes" = "red", "No" = "grey")) +
   coord_equal()

 
 #### save
 outdir <- glue("{plt_dir}/{stamp}")
 dir.create(outdir,showWarnings = F)
 pdf(glue("{outdir}/corr_2.pdf"),width = 6)
 corr + labs(subtitle = "")
 dev.off()
 