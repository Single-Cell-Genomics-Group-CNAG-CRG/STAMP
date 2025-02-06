suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(glue)
  library(here)
  library(SingleCellExperiment)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

# # loading
 dir <- glue("{proj_dir}/data/stamp_7b/raw/stamp_BMP4")
f <- \(.) file.path(dir, paste0("STMPC_1_BMP4ENEC7M_09082024_", .))

y <- readSparseCSV(f("exprMat_file.csv.gz"), transpose=TRUE)
cd <- read.csv(f("metadata_file.csv.gz"))

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
# sce 
sce <- SingleCellExperiment(as, colData=cd, altExps=ae)

library(qs)
qsave(sce, file = glue("{proj_dir}/data/stamp_7b/raw/raw_proc/raw_sce.qs"), nthreads = 8)
