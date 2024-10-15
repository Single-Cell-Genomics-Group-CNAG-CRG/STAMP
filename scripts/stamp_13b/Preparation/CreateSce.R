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
stamp <- "stamp_13b"

# # loading
 dir <- glue("{proj_dir}/data/{stamp}/raw/STAMPC_130K")
 f <- \(.) file.path(dir, paste0("STAMPC_130K_", .))

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
# 
sce <- SingleCellExperiment(as, colData=cd, altExps=ae)

library(qs)
dir.create(glue("{proj_dir}/data/{stamp}/raw/raw_proc"), showWarnings = F,recursive = T)
qsave(sce, file = glue("{proj_dir}/data/{stamp}/raw/raw_proc/raw_sce.qs"), nthreads = 8)