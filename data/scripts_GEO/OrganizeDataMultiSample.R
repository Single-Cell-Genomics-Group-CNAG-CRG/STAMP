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
numb <- "13a"
tech <- "C"
### ### ### ### ### ### ### ### ### ### ### 

# loading
data_dir <- glue("{proj_dir}/data/stamp_13a/raw/raw_proc")
sce <- qread(glue("{data_dir}/layout_sce.qs"))

#dir <- glue("{proj_dir}/data/stamp_{numb}")
#f <- \(.) file.path(dir, paste0("STAMPC_90K", .))

#y <- readSparseCSV(f("exprMat_file.csv.gz"), transpose=TRUE)
#cd <- fread(f("metadata_file.csv.gz"))

# coercion
#y <- as(y[-1, ], "dgCMatrix")
#colnames(y) <- cd$cell

#gs <- rownames(y)
#np <- grep("Negative", gs)
#fc <- grep("SystemControl", gs)

#as <- list(counts=y[-c(np, fc), ])
#ae <- list(
#  negprobes=SingleCellExperiment(list(counts=y[np, ])),
#  falsecode=SingleCellExperiment(list(counts=y[fc, ])))
# 
#sce <- SingleCellExperiment(as, colData=cd, altExps=ae)

#invisible(lapply(names(samples), function(name) {
#  sce$sample[sce$fov %in% samples[[name]]] <<- name
#}))


#sname <- glue("Stamp_{tech}_0{numb}")
#sce <- sce[,! is.na(sce$sample)]
samples_names <- unique(sce$sample)
sname <- glue("Stamp_{tech}_0{numb}")

if(as.numeric(numb) > 10) (sname <- glue("Stamp_{tech}_{numb}")
)

sname <- glue("Stamp_{tech}_{numb}")

samples_names <- samples_names[!is.na(samples_names)]
# Loop through sample names to split & save
lapply(samples_names, \(.){
  
  # Create dirs
  raw_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/CosMx_data/Raw/{sname}_{.}")
  dir.create(raw_dir, showWarnings = F)
  
  proc_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/CosMx_data/Processed/{sname}_{.}")
  dir.create(proc_dir, showWarnings = F)
  
  # Save
  qsave(sce[,sce$sample == .], file = glue("{proc_dir}/{basename({raw_dir})}.qs"), nthreads = 8)
  writeMM(counts(sce[,sce$sample == .]), file = glue("{raw_dir}/{basename({raw_dir})}.mtx"))
  
  library(R.utils)
  
  # Save features
  feature_file <- glue("{raw_dir}/{basename({raw_dir})}_features.tsv")
  fwrite(data.table(features = as.character(rownames(sce[,sce$sample == .]))), quote = FALSE, col.names = FALSE, sep = "\t", file = feature_file)
  gzip(feature_file, overwrite = TRUE)  # Compress to .gz
  
  # Save barcodes
  barcode_file <- glue("{raw_dir}/{basename({raw_dir})}_barcodes.tsv")
  fwrite(data.table(features = as.character(colnames(sce[,sce$sample == .]))), quote = FALSE, col.names = FALSE, sep = "\t", file = barcode_file)
  gzip(barcode_file, overwrite = TRUE)  # Compress to .gz
  
})


