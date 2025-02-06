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

sce <- read10xCounts(glue("{data_dir}/raw/filtered_feature_bc_matrix.h5"),
                     row.names = "symbol",
                     col.names = T,
                     delayed = FALSE,
                     sample.names = "pbmcs_3prime")

counts(sce) <- as(counts(sce),"dgCMatrix")
sub <- downsampleMatrix(counts(sce),prop = 60)

sce$Barcode <- paste0(sce$Sample,"_",colnames(sce))
colnames(sce) <- sce$Barcode

# Save sce object as qs file
dir.create(glue("{data_dir}/raw"), showWarnings = F)

qsave(sce,file = glue("{data_dir}/raw/raw_sce.qs"), nthreads = 8)
