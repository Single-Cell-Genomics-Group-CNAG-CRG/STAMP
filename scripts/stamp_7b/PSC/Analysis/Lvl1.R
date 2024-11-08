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

stamp <- "stamp_7b"
sample <- "iPSCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
sce <- qread(glue("{res_dir}/clustered_sce.qs"))

plotDots(sce, group = "label", features = c("GATA3","IGFBP7","KRT80","KRT19","APOA1","VTN",
                                            "DLL1","SOX2","BMP4","WNT5A","PDGFRA","KDR","TTN","FOXF1",
                                            "SNAI1","WNT3","POU5F1","FGF2"),
         scale = T, center = T) + coord_flip() +
  theme(aspect.ratio = 1/2,
        axis.text.x = element_text(angle = 90, vjust = 1,hjust = 1, color = "black")) 

sce$cluster[sce$label == "b"] <- "pluripotent"
sce$cluster[sce$label == "a"] <- "BMP-induced prog."
sce$cluster[sce$label == "c"] <- "late meso."
sce$cluster[sce$label == "d"] <- "undiff."
sce$cluster[sce$label == "f"] <- "amnion-like"
sce$cluster[sce$label == "e"] <- "ectoderm"

qsave(sce, file = glue("{res_dir}/anno_sce.qs"), nthreads = 8)
