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
})

stamp <- "stamp_18"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/PreProcNew.qs"))
sce <- qread(glue("{res_dir}/qc_sce.qs"))

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

plt_qcmet <- function(var){
  cd$sample <- factor(cd$sample, 
                      levels = names(sort(tapply(cd[[var]], 
                                                 cd$sample, mean, na.rm = TRUE), 
                                          decreasing = TRUE)))
  
  p <- ggplot(cd, aes(y = sample, x = .data[[var]], color = sample, fill = sample)) +
    geom_boxplot(outlier.size = 0.1, alpha = 0.6) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8, color = "black"),
          text = element_text(size = 12, color = "black"),
          aspect.ratio = 2/1,
          legend.position = "right")  +
    labs(fill = "Cell ID", color = "Cell ID")+
    guides(
      fill = guide_legend(ncol = 2),
      color = guide_legend(ncol = 2)
    )
  
  if(var == "sum" | var == "detected") {p <- p + scale_x_log10()}
  if(var == "sum") {p <- p + labs(x = "nCounts")}
  if(var == "detected") {p <- p + labs(x = "nFeature")}
  if(var == "cell_area") {p <- p + labs(x = "Cell area (um2)") + scale_x_continuous(breaks = c(100,300,400))}
  
  return(p)
}

cd <- as.data.frame(colData(sce))
pdf(glue("{outdir}/qcmet_per_sample.pdf"), height = 4, width = 10)
wrap_plots(plt_qcmet("sum") + theme(legend.position = "none", axis.title.y = element_blank()),
           plt_qcmet("detected") + theme(legend.position = "none", axis.title.y = element_blank()),
           plt_qcmet("cell_area") + theme(legend.position = "right",axis.title.y = element_blank()),
           ncol = 3) +
  plot_layout(axis_titles = "collect")
dev.off()