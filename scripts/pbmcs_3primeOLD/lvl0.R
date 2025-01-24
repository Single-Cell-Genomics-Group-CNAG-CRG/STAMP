library(schard)
library(glue)
library(qs)
library(here)
sce <- readH5AD("/Users/emanuelepitino/PhD_Projects/Stamp_V4/data/PBMCs_3prime/raw/INFLAMMATION_main_annotated_LowQFilt_SCGT00_healthy.h5ad")

source(glue("{here()}/scripts/misc/BIN.R")) # bin
source(glue("{here()}/scripts/misc/paths.R")) # bin

proj_dir
dir <- glue("{proj_dir}/data/PBMCs_3prime/raw/INFLAMMATION_main_annotated_LowQFilt_SCGT00_healthy.h5ad")

sce <- h5ad2sce(dir)

sce$level0[sce$Level1 %in% c("T_CD8_NonNaive","T_CD8_Naive")] <- "CD8"
sce$level0[sce$Level1 %in% c("T_CD4_Naive","T_CD4_NonNaive")] <- "CD4"

outdir <- glue("{proj_dir}/data/PBMCs_3prime/processed")
dir.create(outdir,showWarnings = F,recursive = T)
qsave(sce, file = glue("{outdir}/lvl0_sce.qs"))
