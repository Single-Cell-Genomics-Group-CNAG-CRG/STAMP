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
             sce$experiment == "ctrl"]
sub$lvl1 <- as.character(sub$lvl1)
# Proportions
cd <- as.data.frame(colData(sub))

df <- as.data.frame(table(cd$timepoint,cd$lvl1,cd$replicate)) %>%
  group_by(Var1,Var3) %>%
  mutate(pct = round(Freq / sum(Freq),4))

df$Var3 <- as.character(df$Var3)
df$Var3[df$Var3 == "r1"] <- "R.1"
df$Var3[df$Var3 == "r2"] <- "R.2"
df$Var3 <- factor(df$Var3, levels = c("R.1","R.2"))
prop <- ggplot(df, aes(x = Var3, y = pct, fill = Var2)) +
  geom_col() +
#  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(aspect.ratio = 4/1,
        text = element_text(size = 18, color ="black"),
        axis.text = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 18),
        legend.text = element_text(size = 12)) +
  labs(x = "Replicate", y = "Proportion", fill = "Cell type") +
  scale_y_continuous(breaks = c(0.00,0.50,1)) +
  facet_wrap(~Var1) +
  labs(title = unique(sub$experiment))

prop

pdf("/Users/emanuelepitino/Desktop/stamp_13a/Proportions_4h_24h_myeloid_LPS.pdf")
prop
dev.off()
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

plotHeatmap(sub, features = feat,
            order_columns_by  = "sample",
            scale = TRUE,
            center = TRUE)
