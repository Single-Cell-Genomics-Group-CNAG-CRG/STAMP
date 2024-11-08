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

ss <- list(s15=s15, s16=s16, s17=s17, s18=s18)
rm(s15,s16,s17,s18); gc()
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
                "s15","s16","s17","s18")
pal[names(ss)] <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")

# subset shared features
ft <- lapply(ss, rownames)
ft <- Reduce(intersect, ft)
ss <- lapply(ss, \(.) .[ft, ])
rm(s15,s16,s17,s18); gc()
# aggregation
.agg <- \(.) t(rowsum(t(assay(.)), .$sample))
pb <- lapply(ss, .agg)
df <- lapply(pb, \(.) {
  pivot_longer(
    data.frame(gene=rownames(.), .),
    -gene, names_to="cell_line")
}) |> bind_rows(.id="rep")

# correlation
fd <- pivot_wider(df, names_from=c("rep", "cell_line"))
mtx <- cor(as.matrix(select(fd, -gene)), method="pearson")

# plotting
ann <- distinct(df, rep, cell_line)
ann <- data.frame(ann, row.names=rownames(mtx))
col <- list(
  "rep" = pal[names(ss)],
  "cell_line" = pal[1:27])

coln <- colnames(mtx)
ord_cn <- data.frame(stamp_sample = coln)
ord_cn$stamp <- gsub("_.*","", ord_cn$stamp_sample)
ord_cn$sample <- gsub(".*_","", ord_cn$stamp_sample)
ord_cn <- ord_cn %>%
  arrange(ord_cn$sample,ord_cn$stamp)

mtx_2 <- mtx[ord_cn$stamp_sample,ord_cn$stamp_sample]
hm <- pheatmap(mtx_2,
               scale = "none",
               annotation_col = ann,
               annotation_row = ann,
               annotation_colors = col,
               show_colnames = FALSE,
               show_rownames = FALSE,
               border_color = "grey",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               treeheight_row = 0,
               cellwidth = 5,
               cellheight = 5,
               gaps_row = last_indices ,
               gaps_col = last_indices) +
  scale_fill_manual()

hm

stamp <- "high_multi"
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/hm.pdf"),width = 8,height = 7)
hm
dev.off()



# Order points by y-values and get the indices of the top 10 points
top_10_indices <- order(fd[["s18_F"]], decreasing = TRUE)[1:10]

# Plot the points
plot(fd[["s15_F"]], fd[["s18_F"]])
# Add labels from `fd$gene` for the top 10 points
text(fd[["s15_F"]][top_10_indices], fd[["s18_F"]][top_10_indices], 
     labels = fd$gene[top_10_indices], pos = 4, cex = 0.8, col = "blue")
