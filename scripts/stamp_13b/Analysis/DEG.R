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
sce <- sce[,sce$experiment %in% c("LPS","ctrl")]
sce
sce <- logNormCounts(sce)

sce$id <- paste0(sce$lvl1,"_",sce$timepoint,"_",sce$experiment)

# Focus analysis on PDL1+ monocytes (see MDS.R)

sce <- sce[,sce$lvl1 == "PDL1+ mono."]

## Score markers at 4h
sub <- sce[,sce$timepoint == "4h"]
mrk <- scran::scoreMarkers(sub, groups = sub$sample, BPPARAM = bp)

feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(median.logFC.cohen)) %>%
    head(20) %>%
    rownames()
})
feat <- unique(unlist(feat))

plotDots(sub, group = "sample", features = feat, scale = F, center = F) + coord_flip() +
  theme(aspect.ratio = 1/5,
        axis.text.x = element_text(angle = 90, vjust = 1,hjust = 1, color = "black"))



plotExpression(sub, features = "CXCL8", x = "sample", exprs_values = "logcounts")


# Aggregate counts by experiment
c <- counts(sce)  # Take counts matrix
# Identify columns (cells) belonging to each experimental group
lps_cells <- colnames(sce)[sce$experiment == "LPS"]
ctrl_cells <- colnames(sce)[sce$experiment == "ctrl"]
# Rename columns to match their experimental group
colnames(c)[colnames(c) %in% lps_cells] <- "LPS"
colnames(c)[colnames(c) %in% ctrl_cells] <- "ctrl"
# Sum counts for each gene across cells within the same group
c_summed <- t(rowsum(t(c), group = colnames(c), na.rm = TRUE))
# Calculate the number of cells in each group
group_sizes <- c(
  LPS = length(lps_cells),
  ctrl = length(ctrl_cells)
)
# Divide the summed counts by the number of cells to get average counts
c_average <- sweep(c_summed, 2, group_sizes[colnames(c_summed)], FUN = "/")

# Calculate log2 fold change for each gene
log2FC <- log2(c_average[, "LPS"] / c_average[, "ctrl"])
log2FC <- sort(log2FC, decreasing = T)
head(log2FC,30)
