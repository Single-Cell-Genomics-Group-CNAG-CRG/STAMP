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
#sample <- "SKBR3"
### ### ### ### ### ### ### ### ### ### ### 

# Load data

sce <- qread(glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/stamp_{numb}/processed/layout_sce.qs"))
#data_dir <- glue("{proj_dir}/data/stamp_{numb}/{sample}")
#barcodes <- fread(glue("{data_dir}/cell_feature_matrix/barcodes.tsv.gz"), sep = "\t", nThread = 3, header = F)
#features <- fread(glue("{data_dir}/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
#counts <- Matrix::readMM(glue("{data_dir}/cell_feature_matrix/matrix.mtx.gz"))
#md <- fread(glue("{data_dir}/cell_feature_matrix/cells.csv.gz"), sep = ",", nThread = 8, header = T)
counts(sce) <- as(counts(sce),"dgCMatrix")
rownames(colData(sce)) <- sce$cell_id
#colnames(counts) <- barcodes$V1
#rownames(counts) <- features$V2
# split matrix
#gs <- rownames(counts)
#ncp <- grep("NegControlProbe", gs)
#ncc <- grep("NegControlCodeword", gs)
#uc <- grep("UnassignedCodeword", gs)
#dep <- grep("Deprecated",gs)
#int <- grep("Intergenic",gs)
# alternative experiment
#as <- list(counts=counts[-c(ncp, ncc,uc,dep,int), ])
#ae <- list(
#  NegControlProbe=SingleCellExperiment(list(counts=counts[ncp, ])),
#  NegControlCodeword=SingleCellExperiment(list(counts=counts[ncc,])),
#  DeprecatedCodeword = SingleCellExperiment(list(counts = counts[dep,])),
#  IntergenicRegion =  SingleCellExperiment(list(counts = counts[int,])),
#  UnassignedCodeword = SingleCellExperiment(list(counts=counts[uc,])))
# create sce
#sce <- SingleCellExperiment(as, colData=md, altExps=ae)


#sname <- glue("Stamp_{tech}_0{numb}")
#sce$sample <- sample
samples_names <- unique(sce$sample)
sname <- glue("Stamp_{tech}_0{numb}")

if( as.numeric(numb) > 9) (sname <- glue("Stamp_{tech}_{numb}")
)

#samples_names <- samples_names[!is.na(samples_names)]
# Loop through sample names to split & save
lapply(samples_names, \(.){
  
  # Create dirs
  raw_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Xenium_data/Raw/{sname}_{.}")
  dir.create(raw_dir, showWarnings = F)
  
  proc_dir <- glue("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/GEO_submission/Xenium_data/Processed/{sname}_{.}")
  dir.create(proc_dir, showWarnings = F)
  
  # Save
  qsave(sce[,sce$sample == .], file = glue("{proc_dir}/{basename({raw_dir})}.qs"), nthreads = 8)
  writeMM(counts(sce[,sce$sample == .]), file = glue("{raw_dir}/{basename({raw_dir})}.mtx"))
  library(R.utils)
  # Save features
  feature_file <- glue("{raw_dir}/{sname}_{.}_features.tsv")
  fwrite(data.table(features = as.character(rownames(sce[,sce$sample == .]))), quote = FALSE, col.names = FALSE, sep = "\t", file = feature_file)
  gzip(feature_file, overwrite = TRUE)  # Compress to .gz
  
  # Save barcodes
  barcode_file <- glue("{raw_dir}/{sname}_{.}_barcodes.tsv")
  fwrite(data.table(features = as.character(colnames(sce[,sce$sample == .]))), quote = FALSE, col.names = FALSE, sep = "\t", file = barcode_file)
  gzip(barcode_file, overwrite = TRUE)  # Compress to .gz
})