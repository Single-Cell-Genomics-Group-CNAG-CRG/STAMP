suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(qs)
  library(parallel)
  library(scran)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
})

stamp <- "stamp_18"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/PreProcNew.qs"))

# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", use.dimred="PCA", BPPARAM = bp)
# Run Louvain
clusters <- igraph::cluster_louvain(snn.gr, resolution = 1.5)
# Assign to sce
sce$label <- clusters$membership
# Save
qsave(sce, file = glue("{res_dir}/clust_sce.qs"), nthreads = 8)



