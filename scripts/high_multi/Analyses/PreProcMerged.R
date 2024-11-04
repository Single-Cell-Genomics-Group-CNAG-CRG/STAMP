# dependencies
library(ggplot2)
library(dplyr)
library(ggsignif)
library(pheatmap)
library(SingleCellExperiment)
library(glue)
library(here)
library(qs)
library(scater)
library(scuttle)
library(scran)
library(BiocSingular)

stamp <- "high_multi"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
res_dir <- glue("{proj_dir}/data/high_multi/processed")

sce <- qread(glue("{res_dir}/merged.qs"), nthreads = 8)

# subset 25% for each sample
subset_cells <- as.data.frame(colData(sce)) %>%
   mutate(cell_id = rownames(.)) %>%
   group_by(sample) %>% 
   sample_frac(.25) %>%
   pull(cell_id)

sce <- sce[,subset_cells]

counts(sce) <- as(counts(sce), "dgCMatrix") # coerce to dgCMatrix
# LogNormalize 
sce <- logNormCounts(sce, BPPARAM = bp)

sce$s_r_t <- paste0(sce$sample,"_",sce$replicate,"_",sce$tech)

# Feature selection
set.seed(0010101)
dec.var <- modelGeneVar(sce, block = sce$s_r_t, BPPARAM = bp) # model gene var
hvg <- getTopHVGs(dec.var,fdr.threshold = 1) # select hvg on fdr
dec.var$hvg <- "no" # Assign to dec.var column for plot
dec.var$hvg[rownames(dec.var) %in% hvg] <- "yes"
gg_hvg <- plot_hvg(dec.var = dec.var, sub = stamp) # plot
gg_hvg

# PCA
set.seed(101001)
sce <- runPCA(sce, ncomponents=20, BSPARAM=IrlbaParam()) # IRLBA approximation
# run UMAP
set.seed(123)
sce <- runUMAP(sce, dimred="PCA", BPPARAM = bp)

####################################################################################
## PCA plots
####################################################################################
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
pal[28] <- "#A6CEE3"
pal[29] <- "#1F78B4" 
pal[30] <- "#B2DF8A" 
pal[31] <- "#33A02C"
df <- as.data.frame(colData(sce)) # cd
pca <- reducedDim(sce,"PCA") # pc

df <- cbind(df,pca) # merge
df <- df[sample(nrow(df)), ] # shuffle

pca <- \(var){
pal2 <- Polychrome::createPalette(n_distinct(sce[[var]]), c("#FBB4AE", "#B3CDE3", "#CCEBC5"))
names(pal2) <- unique(sce[[var]])


p <- ggplot(df, aes(x = PC1, y = PC2, color = !!sym(var))) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1200) +
  theme_bw() +
  scale_color_manual(values = pal) +
  theme(
    text = element_text(size = 25, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Increase legend dot size
  coord_equal()
if(var == "tech") {p <- p + scale_color_manual(values = pal2)}

return(p)
}


pca_plts <- wrap_plots(
  pca("tech") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  pca("replicate") + theme(legend.position = "none", axis.title.x = element_blank()),
  pca("sample") + theme(legend.position = "none", axis.title.y = element_blank()),
  ncol = 1)

#### save
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/pca_plts.pdf"),width = 6)
pca_plts
dev.off()

qsave(sce, file = glue("{res_dir}/sub_25pct_merged.qs"), nthreads = 8)
