suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(SingleCellExperiment)
  library(glue)
  library(here)
  library(data.table)
  library(qs)
  library(Matrix)
})


dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

### STAMP NUMBER AND TECH ### ### ### ### #
tech <- "Flex"
repl_of <- "Stamp_C_02"
### ### ### ### ### ### ### ### ### ### ### 

# Read the 10X h5 file directly into a sparse matrix
sce <- DropletUtils::read10xCounts("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/flex_cell_lines/sample_filtered_feature_bc_matrix.h5")
rownames(sce) <- rowData(sce)$Symbol
# Add sample metadata if needed
sce$Sample <- NULL
rownames(colData(sce)) <- sce$Barcode
sce$Barcode <- NULL
sce$sample <- "Flex"

counts(sce) <- as(counts(sce),"dgCMatrix")


sname <- glue("Flex_{repl_of}")
savedir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Flex_data")

proc_dir <- glue("{savedir}/Processed/{sname}")
raw_dir <- glue("{savedir}/Raw/{sname}")

dir.create(proc_dir, showWarnings = F)
dir.create(raw_dir, showWarnings = F)

qsave(sce, file = glue("{proc_dir}/{sname}.qs"), nthreads = 8)
writeMM(counts(sce), file = glue("{raw_dir}/{sname}.mtx"))