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

sce$id <-  paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl0 == "myeloid" & 
             sce$timepoint == "4h" &
             sce$experiment %in% c("LPS","ctrl")]

sub <- logNormCounts(sub) # normalize the subset
sub

# Score markers 
mrk <- scran::scoreMarkers(sub, groups = sub$experiment, BPPARAM = bp)
feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(median.logFC.detected)) %>%
    head(15) %>%
    rownames()
})
scoremarkers_feat <- unique(unlist(feat))

sub$sample <- factor(sub$sample, levels = c("ctrl_4h_r1","ctrl_4h_r2",
                                            "LPS_4h_r1","LPS_4h_r2"))

dot_ctrl_lps <- plotDots(sub, group = "sample", features = scoremarkers_feat, scale = TRUE, center = TRUE) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("z-scaled\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.5),
    range = c(1, 5)  # Increase these numbers to make dots bigger
  ) +
  theme_minimal(6) +
  theme(
    aspect.ratio = 1/6,
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 7
    ),
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(),
    legend.key.size = unit(0.5, "lines") 
  ) +
  coord_flip()
dot_ctrl_lps

pdf("/Users/emanuelepitino/Desktop/stamp_13a/dot_myeloid_ctrl_lps.pdf", height = 3)
dot_ctrl_lps
dev.off()
##############################################################################
## 4 vs 24 hours
##############################################################################
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)

sce$id <- paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl0 == "myeloid" & 
             sce$experiment == "LPS"]
sub <- logNormCounts(sub) # normalize the subset
sub

# Score markers 
mrk <- scran::scoreMarkers(sub, groups = sub$sample, BPPARAM = bp)
feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(median.logFC.detected)) %>%
    head(15) %>%
    rownames()
})
scoremarkers_feat <- unique(unlist(feat))

sub$sample <- factor(sub$sample, levels = c("LPS_4h_r1","LPS_4h_r2",
                                            "LPS_24h_r1","LPS_24h_r2"))

dot_myeloid_lps_4h_24h <- plotDots(sub, group = "sample", features = scoremarkers_feat, scale = T, center = T) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("z-scaled\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.5),
    range = c(1, 5)  # Increase these numbers to make dots bigger
  ) +
  theme_minimal(6) +
  theme(
    aspect.ratio = 1/6,
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 7
    ),
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(),
    legend.key.size = unit(0.5, "lines") 
  ) +
  coord_flip()
dot_myeloid_lps_4h_24h
pdf("/Users/emanuelepitino/Desktop/stamp_13a/dot_myeloid_lps_4h_24h.pdf", height = 3)
dot_myeloid_lps_4h_24h
dev.off()








agg <- aggregateAcrossCells(sub, 
                            ids = sub$sample, 
                            use.assay.type = "logcounts",
                            statistics = "mean")

agg$sample <- as.character(agg$sample)

plotHeatmap(agg, features = scoremarkers_feat,
            order_columns_by  = "sample",
            scale = T,
            center = T)

###










sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)

sce$id <-  paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl0 == "myeloid" & 
        #     sce$timepoint == "4h" &
             sce$experiment %in% c("LPS","ctrl")]

sub <- logNormCounts(sub) # normalize the subset
sub$lvl1 <- as.character(sub$lvl1)
sub$sample <- as.character(sub$sample)
agg <- aggregateAcrossCells(sub, 
                            ids = sub$sample, 
                            use.assay.type = "logcounts",
                            statistics = "mean")

agg$sample <- as.character(agg$sample)

plotHeatmap(agg, features = feat,
            order_columns_by  = "sample",
            scale = F,
            center = F) 

pdf("/Users/emanuelepitino/Desktop/stamp_13a/hm_myeloid_log2FC.pdf")
plotHeatmap(agg, features = feat,
            order_columns_by  = "sample",
            scale = TRUE,
            center = TRUE)
dev.off()


# LPS 4h vs 24h
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)
sce$id <-  paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl0 == "myeloid" & 
             sce$experiment == "LPS"]

# Find markers by log2FC of aggregated counts
c <- counts(sub) # Take counts matrix
# Identify columns (cells) belonging to each experimental group
cells_4h <- colnames(sub)[sub$timepoint == "4h"]
cells_24h <- colnames(sub)[sub$timepoint == "24h"]
# Rename columns to match their experimental group
colnames(c)[colnames(c) %in% cells_4h] <- "cells_4h"
colnames(c)[colnames(c) %in% cells_24h] <- "cells_24h"
# Sum counts for each gene across cells within the same group
c_summed <- t(rowsum(t(c), group = colnames(c), na.rm = TRUE))
# Calculate the number of cells in each group
group_sizes <- c(
  cells_4h = length(cells_4h),
  cells_24h = length(cells_24h)
)
# Divide the summed counts by the number of cells to get average counts
c_average <- sweep(c_summed, 2, group_sizes[colnames(c_summed)], FUN = "/")

# Calculate log2 fold change for each gene
log2FC <- log2(c_average[, "cells_4h"] / c_average[, "cells_24h"])
log2FC <- sort(log2FC, decreasing = T)
head(log2FC,30)

feat <- names(log2FC[1:20])
#feat <- unique(unlist(feat))

agg <- aggregateAcrossCells(sub, 
                            ids = sub$sample, 
                            use.assay.type = "logcounts", 
                            statistics = "mean")

agg$sample <- as.character(agg$sample)

plotHeatmap(agg, features = feat,
            order_columns_by  = "sample",
            scale = TRUE,
            center = TRUE)
