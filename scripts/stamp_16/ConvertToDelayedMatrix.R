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
  library(DelayedArray)
})

stamp <- "stamp_16"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/PreProcNew.qs"))
sce

c <- DelayedArray(counts(sce))

sce <- SingleCellExperiment(assay = list(counts = c), colData = colData(sce))

qsave(sce, file = glue("{res_dir}/PreProcNew.qs"), nthreads = 8)
