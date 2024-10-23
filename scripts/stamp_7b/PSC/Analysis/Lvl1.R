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

sce$cluster[sce$label == "d"] <- "amnion-like"
#sce$cluster[sce$label == "a"] <- "mesoderm"
sce$cluster[sce$label == "a"] <- "pluripotent"
sce$cluster[sce$label == "h"] <- "pluripotent"
sce$cluster[sce$label == "g"] <- "BMP-induced prog."
#sce$cluster[sce$label == "d"] <- "mesoderm"
sce$cluster[sce$label == "f"] <- "undiff."
#sce$cluster[sce$label == "f"] <- "mesoderm"
sce$cluster[sce$label == "c"] <- "late meso."
sce$cluster[sce$label == "b"] <- "late meso."
#sce$cluster[sce$label == "g"] <- "mesoderm"
sce$cluster[sce$label == "e"] <- "ectoderm"

qsave(sce, file = glue("{res_dir}/anno_sce.qs"), nthreads = 8)
