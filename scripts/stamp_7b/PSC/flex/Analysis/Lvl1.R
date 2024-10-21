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

sce$cluster[sce$label == "1"] <- "EC1"
sce$cluster[sce$label == "2"] <- "EC1"
sce$cluster[sce$label == "3"] <- "M2"
sce$cluster[sce$label == "4"] <- "EN2"
sce$cluster[sce$label == "5"] <- "EN3"
sce$cluster[sce$label == "6"] <- "EN1"
sce$cluster[sce$label == "7"] <- "P1"
sce$cluster[sce$label == "8"] <- "P2"
sce$cluster[sce$label == "9"] <- "M1"

qsave(sce, file = glue("{outdir}/anno_sce.qs"), nthreads = 8)
