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
tech <- "10X"
repl_of <- "PBMCs"
### ### ### ### ### ### ### ### ### ### ### 
sce <- qread(glue("{proj_dir}/data/PBMCs_5prime/raw/raw_sce.qs"))

# Read the 10X h5 file directly into a sparse matrix
#sce <- DropletUtils::read10xCounts("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/flex_cell_lines/sample_filtered_feature_bc_matrix.h5")
#rownames(sce) <- rowData(sce)$Symbol
# Add sample metadata if needed
#sce$Sample <- NULL
#rownames(colData(sce)) <- sce$Barcode
#sce$Barcode <- NULL
sce$sample <- "5prime"

counts(sce) <- as(counts(sce),"dgCMatrix")


sname <- glue("10X_{repl_of}_{unique(sce$sample)}")
savedir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/10X_data")

proc_dir <- glue("{savedir}/Processed/{sname}")
raw_dir <- glue("{savedir}/Raw/{sname}")

dir.create(proc_dir, showWarnings = F, recursive = T)
dir.create(raw_dir, showWarnings = F, recursive = T)

qsave(sce, file = glue("{proc_dir}/{sname}.qs"), nthreads = 8)
writeMM(counts(sce), file = glue("{raw_dir}/{sname}.mtx"))

library(R.utils)
# Save features
feature_file <- glue("{raw_dir}/{sname}_features.tsv")
fwrite(data.table(features = as.character(rownames(sce))), quote = FALSE, col.names = FALSE, sep = "\t", file = feature_file)
gzip(feature_file, overwrite = TRUE)  # Compress to .gz

# Save barcodes
barcode_file <- glue("{raw_dir}/{sname}_barcodes.tsv")
fwrite(data.table(features = as.character(colnames(sce))), quote = FALSE, col.names = FALSE, sep = "\t", file = barcode_file)
gzip(barcode_file, overwrite = TRUE)  # Compress to .gz
