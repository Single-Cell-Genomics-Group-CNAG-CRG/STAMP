suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(glue)
  library(here)
  library(qs)
  library(SingleCellExperiment)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

stamp <- "stamp_17"
# # loading

sce <- qread( glue("{proj_dir}/data/{stamp}/processed/raw_sce.qs"), nthreads = 8)

layout_plt <- list()
pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

layout <- \(xmin,xmax,region,id1,id2,id3){
ggplot(as.data.frame(colData(sce[,sce$region == region])),
       aes(x = x_centroid, y = y_centroid)) +
  geom_point(shape = 16, size = 0.01) +
  labs(title = region) +
  coord_equal() +
  geom_vline(xintercept = xmin, color = "red") +
  geom_vline(xintercept = xmax, color = "red") 
  
sce$sample[sce$region == region] <<- id2
sce$sample[sce$x_centroid < xmin & sce$region == region] <<- id1
sce$sample[sce$x_centroid > xmax & sce$region == region] <<- id3

plt <- ggplot(as.data.frame(colData(sce[,sce$region == region])),
              aes(x = x_centroid, y = y_centroid)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.001, aes(color = sample)), dpi = 1200) +
  scale_color_manual(values = pal) +
  labs(title = region) +
  coord_equal() + 
  theme(panel.grid = element_blank()) +
  theme_bw() +
  geom_vline(xintercept = xmin, color = "red") +
  geom_vline(xintercept = xmax, color = "red") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = glue("Region: {region}"), color = "Cline ID")
layout_plt[[region]] <<- plt
}

layout(3000,6500,"A","A","H","O")
layout(3000,6500,"B","B","I","P")
layout(3500,6500,"C","C","J","S")
layout(3000,6200,"D","D","K","MX2")
layout(3000,6500,"E","E","L","T")
layout(3000,6200,"F","F","M","U")
layout(3000,6500,"G","G","N","V")

plots_dir <- glue("{plt_dir}/{stamp}")
pdf(glue("{plots_dir}/layout.pdf"), width = 15)
wrap_plots(layout_plt, ncol = 2) +
  plot_layout(axis_titles = "collect")
dev.off()

qsave(sce, file = glue("{proj_dir}/data/{stamp}/processed/layout_sce.qs"), nthreads = 8)
