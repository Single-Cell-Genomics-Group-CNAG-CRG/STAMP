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
dir <- glue("{proj_dir}/data/stamp_{numb}")
f <- \(.) file.path(dir, paste0("STAMPC_90K_", .))

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


#sname <- glue("Stamp_{tech}_0{numb}")
samples_names <- unique(sce$sample)
sname <- glue("Stamp_{tech}_0{numb}")

if(numb > 10) (sname <- glue("Stamp_{tech}_{numb}")
)

samples_names <- samples_names[!is.na(samples_names)]
# Loop through sample names to split & save
lapply(samples_names, \(.){
  
  # Create dirs
  raw_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/CosMx_data/Raw/{sname}_{.}")
  dir.create(raw_dir, showWarnings = F)
  
  proc_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/CosMx_data/Processed/{sname}_{.}")
  dir.create(proc_dir, showWarnings = F)
  
  # Save
  qsave(sce[,sce$sample == .], file = glue("{raw_dir}/{basename({raw_dir})}.mtx"), nthreads = 8)
  writeMM(counts(sce[,sce$sample == .]), file = glue("{proc_dir}/{basename({raw_dir})}.qs"))
})





