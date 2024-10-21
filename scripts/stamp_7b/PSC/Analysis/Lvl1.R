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

sce$cluster[sce$label == "a"] <- "M1"
sce$cluster[sce$label == "b"] <- "P1"
sce$cluster[sce$label == "c"] <- "P2"
sce$cluster[sce$label == "d"] <- "M2"
sce$cluster[sce$label == "e"] <- "EN1"
sce$cluster[sce$label == "f"] <- "M2"
sce$cluster[sce$label == "g"] <- "M2"
sce$cluster[sce$label == "h"] <- "EN2"
sce$cluster[sce$label == "i"] <- "EN3"
sce$cluster[sce$label == "j"] <- "EC1"

qsave(sce, file = glue("{res_dir}/anno_sce.qs"), nthreads = 8)
