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
barcodes <- fread(glue("{data_dir}/filtered_feature_barcode_matrix/A/barcodes.tsv.gz"), sep = "\t", nThread = 3, header = F)
features <- fread(glue("{data_dir}/filtered_feature_barcode_matrix/A/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
counts <- Matrix::readMM(glue("{data_dir}/filtered_feature_barcode_matrix/A/matrix.mtx.gz"))

counts <- as(counts, "dgCMatrix")

colnames(counts) <- barcodes$V1
rownames(counts) <- features$V2

sce <- SingleCellExperiment(assay = list(counts = counts))


# Save sce object as qs file
dir.create(glue("{dir}/raw"), showWarnings = F)

qsave(sce,file = glue("{dir}/raw/raw_sce.qs"), nthreads = 8)
