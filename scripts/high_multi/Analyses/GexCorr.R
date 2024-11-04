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
s16 <- qread(glue("{proj_dir}/data/{stamps[4]}/processed/PreProcNew.qs"))

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


# Function to update counts matrix column names based on sample names
update_counts_colnames <- function(sce) {
  counts_matrix <- assay(sce, "counts")
  sample_names <- setNames(as.character(sce$sample), rownames(colData(sce)))
  colnames(counts_matrix) <- sample_names[colnames(counts_matrix)]
  counts_matrix
}

# Extract count matrix and update colnames with sample names
c17 <- update_counts_colnames(s17)
c18 <- update_counts_colnames(s18)
c15 <- update_counts_colnames(s15)
c16 <- update_counts_colnames(s16)

# function to aggregate and average
aggr_norm <- function(mat){
  # Sum counts across cells grouped by sample name
  summed_counts <- rowsum(t(mat), group = colnames(mat))
  # Compute number of cells per sample
  cell_counts <- table(colnames(mat))
  # Compute average counts per gene for each sample
  average_counts <- sweep(summed_counts, 1, cell_counts, FUN = "/")
  # Transpose back to have genes as rows
  average_counts <- t(average_counts)
  return(average_counts)
}

# aggregate and average
agg_s17 <- as.data.frame(aggr_norm(c17))
agg_s17$gene <- rownames(agg_s17)
agg_s18 <- as.data.frame(aggr_norm(c18))
agg_s18$gene <- rownames(agg_s18)
agg_s15 <- as.data.frame(aggr_norm(c15))
agg_s15$gene <- rownames(agg_s15)
agg_s16 <- as.data.frame(aggr_norm(c16))
agg_s16$gene <- rownames(agg_s16)

# rename columns
colnames(agg_s17)[-which(colnames(agg_s17) == "gene")] <- paste0(colnames(agg_s17)[-which(colnames(agg_s17) == "gene")], "_s17")
colnames(agg_s18)[-which(colnames(agg_s18) == "gene")] <- paste0(colnames(agg_s18)[-which(colnames(agg_s18) == "gene")], "_s18")
colnames(agg_s15)[-which(colnames(agg_s15) == "gene")] <- paste0(colnames(agg_s15)[-which(colnames(agg_s15) == "gene")], "_s15")
colnames(agg_s16)[-which(colnames(agg_s16) == "gene")] <- paste0(colnames(agg_s16)[-which(colnames(agg_s16) == "gene")], "_s16")

df1 <- merge(agg_s17,agg_s18, by = "gene", all = TRUE) # merge
df2 <- merge(agg_s15,agg_s16, by = "gene", all = TRUE) # merge
df <- merge(df1,df2, by = "gene", all = TRUE) # merge
df <- df[,-1] # remove gene column

cor <- cor(df, method = "pearson") # run pearson correlation

anno_col <- data.frame(name = colnames(df))

rownames(anno_col) <- anno_col$name

library(tidyr)
anno_col <- anno_col %>%
  separate(name, into = c("name", "stamp"), sep = "_")

names(anno_col)[names(anno_col) == "name"] <- "Cell line IDs"
names(anno_col)[names(anno_col) == "stamp"] <- "Replicate"
anno_col <- anno_col[, c("Replicate", "Cell line IDs")]
# Create the annotation_colors list
annotation_colors <- list(
  "Replicate" = c("s15" = pal[[28]], "s16" = pal[[29]], "s17" = pal[[30]], "s18" = pal[[31]]),
  "Cell line IDs" = pal[1:27]
)

hm <- pheatmap(
  cor,
  annotation_col = anno_col,
  annotation_row = anno_col,
  scale = "none",
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  treeheight_row = 0,
  treeheight_col = 0
)

stamp <- "high_multi"
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/big_hm.pdf"),width = 8,height = 6)
hm
dev.off()

