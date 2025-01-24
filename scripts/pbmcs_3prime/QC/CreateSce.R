# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(dplyr)
  library(here)
  library(purrr)
  library(glue)
  library(qs)
  library(here)
  library(data.table)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

data_dir <- glue("{proj_dir}/data/PBMCs_3prime")
seu <- readRDS(glue("{data_dir}/raw/DOCTIS_ng5stojs_xsyjzcwr_demuxed_singlets.rds"))

sce <- SingleCellExperiment(assay = list(counts = seu@assays$RNA$counts), colData = seu@meta.data)

# Save sce object as qs file
dir.create(glue("{data_dir}/raw"), showWarnings = F)

qsave(sce,file = glue("{data_dir}/raw/raw_sce.qs"), nthreads = 8)
