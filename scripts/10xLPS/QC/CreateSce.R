# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)
  library(dplyr)
  library(here)
  library(purrr)
  library(glue)
  library(qs)
  library(here)
  library(data.table)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

data_dir <- glue("{proj_dir}/data/10xLPS")

ctrl <- read10xCounts(glue("{data_dir}/raw/ctrl.h5"),
                     row.names = "symbol",
                     col.names = T,
                     delayed = FALSE,
                     sample.names = "ctrl")

lps <- read10xCounts(glue("{data_dir}/raw/lps.h5"),
                      row.names = "symbol",
                      col.names = T,
                      delayed = FALSE,
                      sample.names = "lps")

counts(lps) <- as(counts(lps), "dgCMatrix")
counts(ctrl) <- as(counts(ctrl), "dgCMatrix")

# Save sce object as qs file
dir.create(glue("{data_dir}/raw"), showWarnings = F)
qsave(ctrl,file = glue("{data_dir}/raw/raw_ctrl_sce.qs"), nthreads = 8)
qsave(lps,file = glue("{data_dir}/raw/raw_lps_sce.qs"), nthreads = 8)
