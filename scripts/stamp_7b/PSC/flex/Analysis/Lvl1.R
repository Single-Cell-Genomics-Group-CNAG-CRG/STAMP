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
sample <- "combined"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC/{sample}")
sce <- qread(glue("{outdir}/proc_sce.qs"), nthreads = 8)
sce

plotDots(sce, group = "label", features = c("GATA3","IGFBP7","KRT80","KRT19","APOA1","VTN",
                                            "DLL1","SOX2","BMP4","WNT5A","PDGFRA","KDR","TTN","FOXF1",
                                            "SNAI1","WNT3","POU5F1","FGF2"),
         scale = T, center = T) + coord_flip() +
  theme(aspect.ratio = 1/2,
        axis.text.x = element_text(angle = 90, vjust = 1,hjust = 1, color = "black"))

sce$cluster[sce$label == "1"] <- "ectoderm"
sce$cluster[sce$label == "2"] <- "ectoderm"
sce$cluster[sce$label == "3"] <- "ectoderm"
sce$cluster[sce$label == "4"] <- "ectoderm"
sce$cluster[sce$label == "5"] <- "ectoderm"

sce$cluster[sce$label == "7"] <- "pluripotent"
sce$cluster[sce$label == "8"] <- "pluripotent"
sce$cluster[sce$label == "9"] <- "pluripotent"

sce$cluster[sce$label == "10"] <- "amnion-like"
sce$cluster[sce$label == "11"] <- "pluripotent"
sce$cluster[sce$label == "12"] <- "pluripotent"

qsave(sce, file = glue("{outdir}/anno_sce.qs"), nthreads = 8)
