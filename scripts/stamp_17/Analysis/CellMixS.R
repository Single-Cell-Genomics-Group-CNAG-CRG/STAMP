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
  library(reshape2)
})

stamp <- "stamp_17"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/clust_sce.qs"))

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")




# Load required packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(cowplot)
  library(limma)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scater)
  library(CellMixS)
})

# Load sim_list example data
sim_list <- readRDS(system.file(file.path("extdata", "sim50.rds"), 
                                package = "CellMixS"))
names(sim_list)
#> [1] "batch0"  "batch20" "batch50"

sce50 <- sim_list[["batch50"]]

sub <- sce[,sample(colnames(sce),100000)]
sce50 <- cms(sub, k = 30, group = "sample", res_name = "unaligned", 
             n_dim = 3, cell_min = 4)
visHist(sce50)
