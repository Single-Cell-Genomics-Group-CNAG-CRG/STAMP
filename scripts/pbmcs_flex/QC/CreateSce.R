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

data_dir <- glue("{proj_dir}/data/flex_PBMCs")

sce <- read10xCounts(glue("{data_dir}/raw/sample_filtered_feature_bc_matrix.h5"),
                     row.names = "symbol",
                     col.names = T,
                     delayed = FALSE,
                     sample.names = "flex")

counts(sce) <- as(counts(sce), "dgCMatrix")

sce <- sce[,sample(colnames(sce),10000)] # subset to keep 10k cells
# Save sce object as qs file
dir.create(glue("{data_dir}/raw"), showWarnings = F)
qsave(sce,file = glue("{data_dir}/raw/raw_sce.qs"), nthreads = 8)
