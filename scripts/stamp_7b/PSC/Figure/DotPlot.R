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
#cosmx <- qread(glue("{res_dir}/clustered_sce.qs"))
cosmx$tech <- "cosmx"

cosmx$cluster <- factor(cosmx$cluster,
                        levels = c("amnion-like",
                                   "BMP-induced prog.",
                                   "ectoderm",
                                   "late meso.",
                                   "pluripotent",
                                   "undiff."))

feat <-c("GATA3","IGFBP7","KRT80","KRT19",
         "DLL1","BMP4",
         "SOX2","NRG1",
         "WNT5A","PDGFRA","KDR","TTN","FOXF1",
         "SOX2","NRG1",
         "POU5F1","FGF2","SNAI1",
         "EOMES",
         "WNT3")

feat <- factor(feat, levels = c("GATA3","IGFBP7","KRT80","KRT19",
                                "DLL1","SOX2","NRG1",
                                "BMP4","WNT5A","PDGFRA","KDR",
                                "TTN","FOXF1","SNAI1",
                                "EOMES",
                                "WNT3","POU5F1","FGF2"))

gg_dot <- plotDots(cosmx,feat, group = "cluster",
         scale = T, center = T) +
  scale_color_gradient2("z-scaled\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.5),
    range = c(2, 6)  # Increase these numbers to make dots bigger
  ) +
  theme_minimal(6) +
  theme(
    aspect.ratio = 1/2,
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 12
    ),
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(),
    legend.key.size = unit(0.5, "lines") 
  ) +
  coord_flip()

gg_dot

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC/dot.pdf")
gg_dot
dev.off()
