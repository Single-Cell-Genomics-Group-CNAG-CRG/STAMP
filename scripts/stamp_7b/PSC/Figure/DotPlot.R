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
sce <- qread(glue("{res_dir}/PreProcNew.qs"))


feat <- c("FGF2","POU5F1","WNT3","KDR",
          "GATA3","KRT19","FOXF1","BMP4","IGFBP7","TTN","EOMES",
          "SNAI1","DLL1","SOX2","NRG1")
feat <- factor(feat, levels = feat)

gg_dot <- plotDots(sce,feat, group = "sample",
         scale = T, center = T) +
  scale_color_gradient2("z-scaled\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.5),
    range = c(2, 8)  # Increase these numbers to make dots bigger
  ) +
  theme_minimal(6) +
  theme(
    aspect.ratio = 1/3,
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

pdf("/Users/emanuelepitino/Desktop/stamp_7b_PSC_V2/dot.pdf", width = 6, height = 3)
gg_dot
dev.off()
