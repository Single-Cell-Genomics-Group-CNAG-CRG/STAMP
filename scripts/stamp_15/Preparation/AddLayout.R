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

stamp <- "stamp_13a"
# # loading

sce <- qread( glue("{proj_dir}/data/{stamp}/raw/raw_proc/raw_sce.qs"), nthreads = 8)


# Add samples information from layout and artifacts from AtoMx
samples <- list(
  A = c(1:30),
  B = c(31:60),
  C = c(61:90),
  D = c(91:120),
  E = c(121:150),
  F = c(151:180),
  
  G = c(181:210),
  H = c(211:240),
  I = c(241:270),
  J = c(271:300),
  K = c(301:330),
  L = c(331:360),
  
  M = c(361:390),
  N = c(391:420),
  O = c(421:450),
  P = c(451:480),
  Q = c(481:510),
  R = c(511:540),
  
  S = c(541:570),
  T = c(571:600),
  U = c(601:630),
  V = c(631:660),
  W = c(661:690),
  MX1 = c(691:720)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

qsave(sce, file = glue("{proj_dir}/data/{stamp}/raw/raw_proc/layout_sce.qs"), nthreads = 8)
