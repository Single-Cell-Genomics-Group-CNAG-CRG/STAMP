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

stamp <- "stamp_16"
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
  
  G = c(181:205),
  H = c(206:235),
  I = c(236:265),
  J = c(266:290),
  K = c(291:315),
  L = c(316:340),
  
  M = c(341:370),
  N = c(371:395),
  O = c(396:420),
  P = c(421:450),
  Q = c(451:475),
  R = c(476:500),
  
  S = c(501:525),
  T = c(526:555),
  U = c(556:585),
  V = c(586:615),
  X = c(616:640),
  Y = c(641:665)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

qsave(sce, file = glue("{proj_dir}/data/{stamp}/raw/raw_proc/layout_sce.qs"), nthreads = 8)
