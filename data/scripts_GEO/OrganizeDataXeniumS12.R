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
numb <- "12"
tech <- "X"
### ### ### ### ### ### ### ### ### ### ### 

sub <- "PBMCs"
# Load data
data_dir <- glue("{proj_dir}/data/stamp_{numb}/{sub}")
barcodes <- fread(glue("{data_dir}/cell_feature_matrix/barcodes.tsv.gz"), sep = "\t", nThread = 3, header = F)
features <- fread(glue("{data_dir}/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
counts <- Matrix::readMM(glue("{data_dir}/cell_feature_matrix/matrix.mtx.gz"))
md <- data.table::fread(glue("{data_dir}/cell_feature_matrix/cells.csv.gz"))

counts <- as(counts, "dgCMatrix")

colnames(counts) <- barcodes$V1
rownames(counts) <- features$V2

gs <- rownames(counts)
ncp <- grep("NegControlProbe", gs)
ncc <- grep("NegControlCodeword", gs)
uc <- grep("UnassignedCodeword", gs)

as <- list(counts=counts[-c(ncp, ncc,uc), ])
ae <- list(
  NegControlProbe=SingleCellExperiment(list(counts=counts[ncp, ])),
  NegControlCodeword=SingleCellExperiment(list(counts=counts[ncc,])),
  UnassignedCodeword = SingleCellExperiment(list(counts=counts[uc,])))

sce <- SingleCellExperiment(as, colData = md, altExps=ae)


sname <- glue("Stamp_{tech}_0{numb}_{sub}")
if(as.numeric(numb) > 9) (sname <- glue("Stamp_{tech}_{numb}_{sub}")
)
sce
savedir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Xenium_data")

proc_dir <- glue("{savedir}/Processed/{sname}")
raw_dir <- glue("{savedir}/Raw/{sname}")

dir.create(proc_dir, showWarnings = F)
dir.create(raw_dir, showWarnings = F)

qsave(sce, file = glue("{proc_dir}/{sname}.qs"), nthreads = 8)
writeMM(counts(sce), file = glue("{raw_dir}/{sname}.mtx"))