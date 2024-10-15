suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(glue)
  library(here)
  library(qs)
  library(SingleCellExperiment)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

stamp <- "stamp_13b"
# # loading

sce <- qread( glue("{proj_dir}/data/{stamp}/raw/raw_proc/raw_sce.qs"), nthreads = 8)


# Add samples information from layout and artifacts from AtoMx
samples <- list(
  ctrl_24h_r1 = c(226:250),
  ctrl_24h_r2 = c(151:175),
  ctrl_4h_r1 = c(76:100),
  ctrl_4h_r2 = c(1:25),
  
  LPS_24h_r1 = c(251:275),
  LPS_24h_r2 = c(176:200),
  LPS_4h_r1 = c(101:125),
  LPS_4h_r2 = c(26:50),
  
  aCD3aCD28_24h_r1 = c(276:300),
  aCD3aCD28_24h_r2 = c(201:225),
  aCD3aCD28_4h_r1 = c(126:150),
  aCD3aCD28_4h_r2 = c(51:75)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

qsave(sce, file = glue("{proj_dir}/data/{stamp}/raw/raw_proc/layout_sce.qs"), nthreads = 8)
