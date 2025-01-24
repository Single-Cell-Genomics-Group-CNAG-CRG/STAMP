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
repl_of <- "Stamp_C_07b"
### ### ### ### ### ### ### ### ### ### ### 



base_dir <- "/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/flex_stamp_7b"

# Find the folder that ends with the value in `sample`
dir <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

sce_ll <- lapply(dir, \(.){

  mtx <- Matrix::readMM(glue("{.}/matrix.mtx.gz"))
  mtx <- as(mtx, "dgCMatrix") # convert to dgCMatrix
  bc <- fread(glue("{.}/barcodes.tsv.gz"), header = F)
  feat <- fread(glue("{.}/features.tsv.gz"), header = F)
  rownames(mtx) <- feat$V2
  colnames(mtx) <- bc$V1
  
  sce <- SingleCellExperiment(assays = list(counts = mtx), rowData = feat)
  
  # Define the file path
  file_path <- .
  
  # Extract the specific string after the last underscore
  samp <- sub(".*_(\\w+)$", "\\1", file_path)
  sce$sample <- samp
  return(sce)
})

names(sce_ll) <- sapply(sce_ll, \(.) unique(.$sample))

samples <- names(sce_ll)

for(samp in samples){
  
sname <- glue("Flex_{repl_of}_{samp}")
savedir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Flex_data")

proc_dir <- glue("{savedir}/Processed/{sname}")
raw_dir <- glue("{savedir}/Raw/{sname}")

dir.create(proc_dir, showWarnings = F)
dir.create(raw_dir, showWarnings = F)

qsave(sce_ll[[samp]], file = glue("{proc_dir}/{sname}.qs"), nthreads = 8)
writeMM(counts(sce_ll[[samp]]), file = glue("{raw_dir}/{sname}.mtx"))
}
