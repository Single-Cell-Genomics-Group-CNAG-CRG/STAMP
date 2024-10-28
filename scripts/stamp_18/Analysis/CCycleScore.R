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
sce <- qread(glue("{res_dir}/PreProcNew.qs"))
sce <- qread(glue("{res_dir}/qc_sce.qs"))
rownames(colData(sce)) <- sce$cell_id

pal <- Polychrome::createPalette(21, c("#8DD3C7","#BEBADA","#FB8072"))
names(pal) <- c("A","O","H",
                "B","I","P",
                "C","J","S",
                "D","K","MX2",
                "E","L","T",
                "F","M","U",
                "G","N","V")

#sub <- sce[,sample(colnames(sce),10000)]
sub <- sce

pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                package="scran"))

assignments <- cyclone(sub, pairs, BPPARAM = bp, gene.names = rowData(sub)$ID)


sub$phase <- assignments$phases

cd <- as.data.frame(colData(sub))

# Assuming 'cd' is your data frame with 'sample' and 'phase' columns
# Compute the counts table
tab <- table(cd$sample, cd$phase)

# Calculate the proportions
prop_tab <- prop.table(tab, margin = 1)

# Convert the table to a data frame
df <- as.data.frame(as.table(prop_tab))
colnames(df) <- c("Sample", "Phase", "Proportion")

# Reshape data to wide format for clustering
df_wide <- dcast(df, Sample ~ Phase, value.var = "Proportion")

# Set row names and remove the 'Sample' column
rownames(df_wide) <- df_wide$Sample
df_wide$Sample <- NULL

# Compute the distance matrix (Euclidean distance)
dist_matrix <- dist(df_wide)

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Get the order of samples from the clustering result
sample_order <- hc$labels[hc$order]

# Reorder the levels of 'Sample' according to the clustering
df$Sample <- factor(df$Sample, levels = sample_order)

# Ensure 'Phase' is an ordered factor if needed
df$Phase <- factor(df$Phase, levels = c("G2M", "S", "G1"))

# Plot using ggplot2
ccycle <- ggplot(df, aes(y = Sample, x = Proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    text = element_text(size = 18, color = "black"),
    aspect.ratio = 2 / 1
  ) +
  labs(y = "Sample", x = "Proportion", fill = "Cell cycle \nPhase")

outdir <- glue("{plt_dir}/{stamp}/QC")
pdf(glue("{outdir}/ccycle.pdf"), height = 4, width = 10)
ccycle
dev.off()
  