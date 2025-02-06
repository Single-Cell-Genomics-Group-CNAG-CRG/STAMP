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
numb <- "1"
tech <- "C"
### ### ### ### ### ### ### ### ### ### ### 

# loading
dir <- glue("{proj_dir}/data/stamp_{numb}")
f <- \(.) file.path(dir, paste0("PBMCs_", .))

y <- readSparseCSV(f("exprMat_file.csv.gz"), transpose=TRUE)
cd <- fread(f("metadata_file.csv.gz"))

# coercion
y <- as(y[-1, ], "dgCMatrix")
colnames(y) <- cd$cell

gs <- rownames(y)
np <- grep("Negative", gs)
fc <- grep("SystemControl", gs)

as <- list(counts=y[-c(np, fc), ])
ae <- list(
  negprobes=SingleCellExperiment(list(counts=y[np, ])),
  falsecode=SingleCellExperiment(list(counts=y[fc, ])))
# 
sce <- SingleCellExperiment(as, colData=cd, altExps=ae)



sname <- glue("Stamp_{tech}_0{numb}")
savedir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/CosMx_data")

dir.create(glue("{savedir}/Processed/{sname}"))
dir.create(glue("{savedir}/Raw/{sname}"))
qsave(sce, file = glue("{savedir}/Processed/{sname}/{sname}.qs"), nthreads = 8) # save sce obj
writeMM(counts(sce), file = glue("{savedir}/Raw/{sname}/{sname}.mtx")) # save matrix

library(R.utils)

# Save features
feature_file <- glue("{savedir}/Raw/{sname}/features.tsv")
fwrite(data.table(features = as.character(rownames(sce))), quote = FALSE, col.names = FALSE, sep = "\t", file = feature_file)

# Save barcodes
barcode_file <- glue("{savedir}/Raw/{sname}/barcodes.tsv")
fwrite(data.table(features = as.character(colnames(sce))), quote = FALSE, col.names = FALSE, sep = "\t", file = barcode_file)
