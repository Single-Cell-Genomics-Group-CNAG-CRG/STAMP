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
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

dir <- glue("{proj_dir}/data/PBMCs_5prime")


sce <- read10xCounts(glue("{dir}/raw/sample_filtered_feature_bc_matrix.h5"),
                       row.names = "symbol",
                       col.names = T,
                       delayed = FALSE,
                       sample.names = "5prime")

counts(sce) <- as(counts(sce),"dgCMatrix")

sce$Barcode <- paste0(sce$Sample,"_",colnames(sce))
colnames(sce) <- sce$Barcode
  

# Save sce object as qs file
dir.create(glue("{dir}/raw"), showWarnings = F)

qsave(sce,file = glue("{dir}/raw/raw_sce.qs"), nthreads = 8)