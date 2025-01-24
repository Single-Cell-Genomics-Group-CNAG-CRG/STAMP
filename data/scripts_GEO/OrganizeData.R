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
numb <- "2"
tech <- "C"
### ### ### ### ### ### ### ### ### ### ### 

# loading
dir <- glue("{proj_dir}/data/stamp_{numb}")
f <- \(.) file.path(dir, paste0("SML_Square_", .))

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

qsave(sce, file = glue("{savedir}/Processed/{sname}/{sname}.qs"), nthreads = 8)
writeMM(counts(sce), file = glue("{savedir}/Raw/{sname}/{sname}.mtx"))




