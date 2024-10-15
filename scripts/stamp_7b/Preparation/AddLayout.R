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

# # loading
dir <- glue("{proj_dir}/data/stamp_7b/raw/raw_proc")

sce <- qread( glue("{proj_dir}/data/stamp_7b/raw/raw_proc/raw_sce.qs"), nthreads = 8)


# Add samples information from layout and artifacts from AtoMx
samples <- list(
  iESC_0h = c(1:64),
  iESC_6h = c(129:189),
  iESC_12h = c(251:310),
  iESC_24h = c(546:600),
  iESC_48h = c(65:128),
  iESC_72h = c(190:250),
  iESC_96h = c(311:372),
  iESC_120h = c(601:654),
  ectoderm = c(490:545),
  mesoderm = c(434:489),
  endoderm = c(373:433),
  iPSC_parental = c(655:708)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

qsave(sce, file = glue("{proj_dir}/data/stamp_7b/raw/raw_proc/layout_sce.qs"), nthreads = 8)
