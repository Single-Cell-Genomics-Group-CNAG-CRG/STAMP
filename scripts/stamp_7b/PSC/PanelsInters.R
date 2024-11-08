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

# Load CosMx data
res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
cosmx <- qread(glue("{res_dir}/qc_sce.qs"))

cosmx <- cosmx[,cosmx$sample != "endoderm"]
cosmx

# Load Flex data
base_dir <- glue("{proj_dir}/data/{stamp}/processed/flex/iPSC_parental")

flex <- qread(glue("{base_dir}/qc_flex_iPSC_parental.qs"), nthreads = 8)
flex

feat <- intersect(rownames(cosmx),rownames(flex))

qsave(feat, file = glue("{proj_dir}/data/{stamp}/processed/panel_int.qs"))


