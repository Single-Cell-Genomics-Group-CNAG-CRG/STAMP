# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(pheatmap)
library(SingleCellExperiment)
library(glue)
library(here)
library(qs)

stamps <- c("stamp_17","stamp_18","stamp_15","stamp_16")
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# load data
s17 <- qread(glue("{proj_dir}/data/{stamps[1]}/processed/clust_sce.qs"))
s18 <- qread(glue("{proj_dir}/data/{stamps[2]}/processed/clust_sce.qs"))
s15 <- qread(glue("{proj_dir}/data/{stamps[3]}/processed/PreProcNew.qs"))

s16 <- s15 # mock, remove once s16 is ready
set.seed(123)
pal <- Polychrome::createPalette(31, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V",
                "W","X","Y","MX1","Q","R",
                "stamp_15","stamp_16","stamp_17","stamp_18")
pal[28] <- "#A6CEE3"
pal[29] <- "#1F78B4" 
pal[30] <- "#B2DF8A" 
pal[31] <- "#33A02C"

ft <- intersect(rownames(s15),rownames(s18))
s15 <- s15[ft,]
s16 <- s16[ft,]
s17 <- s17[ft,]
s18 <- s18[ft,]

sub <- \(sce){sce[,sample(colnames(sce),10000)]}

fun_hvg <- \(sce) {
  sce <- logNormCounts(sce)
  dec.var <- modelGeneVar(sce, block = sce$sample, BPPARAM = bp)
  hvg <- getTopHVGs(dec.var, fdr.threshold = 1)
  return(hvg)
}

s17 <- sub(s17)
s15 <- sub(s15)

hvg_s17 <- fun_hvg(s17)
hvg_s15 <- fun_hvg(s15)


x <- list(
  CosMx = hvg_s15, 
  Xenium = hvg_s17
)

ggvenn(
  x, 
  fill_color = c("#A6CEE3","#B2DF8A"),
  stroke_size = 0.5, set_name_size = 4
)


length(intersect(hvg_s15,hvg_s17)) / n_distinct(c(hvg_s15,hvg_s17))