
library(zellkonverter)
sce <- readH5AD("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/PBMCs_3prime/raw/INFLAMMATION_main_annotated_LowQFilt_SCGT00_healthy.h5ad")

dir <- "/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/PBMCs_3prime/raw/INFLAMMATION_main_annotated_LowQFilt_SCGT00_healthy.h5ad"
library(schard)
sce <- h5ad2sce(dir)

sce$level0[sce$Level1 %in% c("T_CD8_NonNaive","T_CD8_Naive")] <- "CD8"
sce$level0[sce$Level1 %in% c("T_CD4_Naive","T_CD4_NonNaive")] <- "CD4"
