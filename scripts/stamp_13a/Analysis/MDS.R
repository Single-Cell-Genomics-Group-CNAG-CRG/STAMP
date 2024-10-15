suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(here)
  library(scuttle)
  library(glue)
  library(qs)
})

# Data loading
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_13a"
res_dir <- glue("{proj_dir}/data/{stamp}/processed")
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)
sce <- sce[,sce$experiment %in% c("aCD3aCD28","ctrl")]
sce
sce <- logNormCounts(sce)

sce$id <- paste0(sce$replicate,"_",sce$lvl1,"_",sce$timepoint,"_",sce$experiment)

#agg <- aggregateAcrossCells(sce, ids = sce$id, use.assay.type = "logcounts", statistics = "mean")
agg <- aggregateAcrossCells(sce, ids = sce$sample, use.assay.type = "logcounts")
agg <- runMDS(agg)

mds_plt <- list()
for(sub in unique(agg$lvl1)){
  
  mds_plt[[sub]]<- plotMDS(agg[,agg$lvl1 == sub], 
                           size_by = "replicate", color_by = "experiment", shape_by = "timepoint", point_size =5) +
    labs(subtitle = sub) +
    theme(plot.subtitle = element_text(size = 15, color = "black"))
}


plotMDS(agg, color_by = "experiment", shape_by = "replicate", size_by = "timepoint")

wrap_plots(mds_plt) + plot_layout(guides = "collect")wrap_plots(mds_plt) + plot_lvl1layout(guides = "collect")
