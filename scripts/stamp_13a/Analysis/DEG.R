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
    arrange(desc(mean.logFC.cohen)) %>%
    head(15) %>%
    rownames()
})

scoremarkers_feat <- unique(unlist(feat))

sub$sample <- factor(sub$sample, levels = c("ctrl_4h_r1","ctrl_4h_r2",
                                            "LPS_4h_r1","LPS_4h_r2"))

dot_ctrl_lps <- plotDots(sub, group = "sample", features = scoremarkers_feat,
                         scale = F, center = T) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("centered\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
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
    arrange(desc(mean.logFC.cohen)) %>%
    head(15) %>%
    rownames()
})
scoremarkers_feat <- unique(unlist(feat))

sub$sample <- factor(sub$sample, levels = c("LPS_4h_r1","LPS_4h_r2",
                                            "LPS_24h_r1","LPS_24h_r2"))

dot_myeloid_lps_4h_24h <- plotDots(sub, group = "sample", features = scoremarkers_feat,
                                   scale = F, center = T) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("centered\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
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
###


sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)

sce$id <- paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl0 == "myeloid" & 
             sce$experiment %in% c("LPS","ctrl")]
sub <- logNormCounts(sub) # normalize the subset
sub

# Score markers 
mrk <- scran::scoreMarkers(sub, groups = sub$sample, BPPARAM = bp)

feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(mean.logFC.cohen)) %>%
    head(15) %>%
    rownames()
})
scoremarkers_feat <- unique(unlist(feat))

feats <- c("CD80","CD86","IL1B","TLR1","TLR2","TLR3","TLR4","TLR5",
           "TLR7","TLR8")

plotDots(sub, group = "sample", features = feats,
         scale = F, center = T) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("centered\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
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


###
sce <- qread(glue("{res_dir}/anno_sce_P1.qs"), nthreads = 8)

sce$id <- paste0(sce$lvl1,"_",sce$replicate,"_",sce$timepoint,"_",sce$experiment)
sce

# subset to keep only myeloid compartment
sub <- sce[,sce$lvl1 == "PDL1+ mono." & 
             sce$experiment %in% c("ctrl","LPS")]

sub <- logNormCounts(sub) # normalize the subset
sub

# Score markers 
mrk <- scran::scoreMarkers(sub, groups = sub$sample, BPPARAM = bp)

feats <- c("TLR1","TLR2","TLR3","TLR4","TLR5",
           "TLR7","TLR8","IL1B","IL6","IL10","TNF","TGFBI","NFKB","
           CD80","CD86","CXCL3","CCL3","CCL2","CXCL8","S100A8","S100A9",
           "CCL20","IDO1")

feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(mean.logFC.cohen)) %>%
    head(15) %>%
    rownames()
})
scoremarkers_feat <- unique(unlist(feat))

plotDots(sub, group = "sample", features = intersect(feats, rownames(sub)),
         scale = F, center = T) +
  theme(aspect.ratio = 3/1) +
  scale_color_gradient2("centered\nmean expr.", low = "blue4", mid = "grey90", high = "red4") +
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
