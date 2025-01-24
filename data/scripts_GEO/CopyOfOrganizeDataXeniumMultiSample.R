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
numb <- "18"
tech <- "X"
### ### ### ### ### ### ### ### ### ### ### 

# Load data

sce <- qread("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/stamp_18/processed/layout_sce.qs")

counts(sce) <- as(counts(sce),"dgCMatrix")
#sname <- glue("Stamp_{tech}_0{numb}")
samples_names <- unique(sce$sample)
sname <- glue("Stamp_{tech}_0{numb}")

if( as.numeric(numb) > 9) (sname <- glue("Stamp_{tech}_{numb}")
)

samples_names <- samples_names[!is.na(samples_names)]
# Loop through sample names to split & save
lapply(samples_names, \(.){
  
  # Create dirs
  raw_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Xenium_data/Raw/{sname}_{.}")
  dir.create(raw_dir, showWarnings = F)
  
  proc_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Xenium_data/Processed/{sname}_{.}")
  dir.create(proc_dir, showWarnings = F)
  
  # Save
  qsave(sce[,sce$sample == .], file = glue("{raw_dir}/{basename({raw_dir})}.mtx"), nthreads = 8)
  writeMM(counts(sce[,sce$sample == .]), file = glue("{proc_dir}/{basename({raw_dir})}.qs"))
})