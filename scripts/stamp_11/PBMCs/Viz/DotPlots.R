# dependencies
library(qs)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(glue)
library(here)

dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_11"
sub <- "PBMCs"


# B
################################################################################
################################################################################
lin <- "B"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/{lin}/lvl2_sce.qs"), nthreads = 8)

b <- c("CD79A","CD79B","CD19","IGHM","CD27","CD24","CD38","ITGAX","SELL","CCR7","TCL1A")
b2 <- c("BCL6", "CD19", "CD79A", "CD79B","CD27","CD38", "CD40", "CXCR5", "MS4A1", "MYC")

fs <- intersect(unique(c(b,b2)), rownames(sce))

sce$lvl2 <- factor(sce$lvl2, levels = c("transitional B","B activated - plasma cells", "B naive", "B mature", "B memory"))
fs <- c("CD19","TCL1A","CD79A","CD79B","MYC","CD38","CXCR5","MZB1","SELL","CD40","CCR7","BCL6","CD27","ITGAX","MS4A1",
        "MZB1")

fs <- intersect(fs, rownames(sce))

dot_b <- plotDots(sce, features = fs, group = "lvl2", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  coord_flip() + 
  theme_minimal(6) + 
  theme(
    aspect.ratio = length(unique(sce$lvl2))/length(fs),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text = element_text(size = 14, color = "black")
  ) + 
  scale_size_continuous(range = c(3, 7), breaks = c(0.25, 0.75)) 

dot_b

b <- sce
# Myeloid
################################################################################
################################################################################
lin <- "myeloid"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/{lin}/lvl2_sce.qs"), nthreads = 8)

fs <- c("IL1B","FCGR3A","CD14","CD1C","CCR2","CD68","ITGAM","CD34","FLT3","CX3CR1","CLEC9A","CCR5",
          "S100A8","S100A9","LYZ","ITGAX")

fs <- intersect(unique(fs), rownames(sce))



sce$lvl2 <- factor(sce$lvl2, levels = c("class. mono.", "infl. mono.", "int. mono.","non-class. mono.",
                                        "DC4","DC1/DC2/DC5","pDC"))


fs <- c("CD14","ITGAM","IL1B","CD68","FCGR3A","CX3CR1","ITGAX","CD1C","CD34","CLEC9A","FLT3","CCR5","CCR2")

dot_myelo <- plotDots(sce, features = fs, group = "lvl2", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  coord_flip() + 
  theme_minimal(6) + 
  theme(
    aspect.ratio = length(unique(sce$lvl2))/length(fs),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text = element_text(size = 14, color = "black")
  ) + 
  scale_size_continuous(range = c(3, 10), breaks = c(0.25, 0.75)) 

dot_myelo

myelo <- sce
# CD8
################################################################################
################################################################################
lin <- "CD8"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/{lin}/lvl2_sce.qs"), nthreads = 8)

t <- c("CCR5","GZMK", "EOMES","LAG3","ZNF683","GZMH","KLRD1","PDCD1","FASLG","CCR6",
       "CCR4","CCL20","RORC","CCR7","KLRK1","REG4","CD38","FOXP3","IL7R","IL12RB2",
       "CD40LG","SELL","LEF1","KLRB1")


cd8 <- c("TCF7", "KLRB1", "KLRC1", "KLRD1", "KLRF1", "KLRG1", "CX3CR1", "ENTPD1", "IFI16", "IFI35", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "GZMA", "GZMB", "GZMH", "GZMK", "HAVCR2")


fs <- c("CD38","TNFRSF9","HAVCR2","PDCD1","CCR5", "GZMK", "EOMES", "LAG3","KLRG1","IFI35","IFIH1",
        "CCR6","CCL20","RORC","IL12RB2","CD40LG",
        "KLRB1","KLRC1","ENTPD1","IFI44L","CCR4","FOXP3",
        "SELL","LEF1","TCF7","CCR7","IFITM1","CX3CR1","GZMA","ZNF683", "GZMH","GZMB", "KLRD1",
        "FASLG")

sce$lvl2 <- factor(sce$lvl2, levels = c("CD8 T act.","CD8 T inf.resp.","CD8 TCM","CD8 TN","CD8 T eff.","CD8 TEM",
                                        "CD8 TPM","CD8 T cyto."))



  
dot_cd8 <- plotDots(sce, features = fs, group = "lvl2", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  coord_flip() + 
  theme_minimal(6) + 
  theme(
    aspect.ratio = length(unique(sce$lvl2))/length(fs),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text = element_text(size = 8, color = "black")
  ) + 
  scale_size_continuous(range = c(2, 4), breaks = c(0.25, 0.75)) 

dot_cd8

cd8 <- sce
# CD4
################################################################################
################################################################################
lin <- "CD4"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/{lin}/lvl2_sce.qs"), nthreads = 8)

fs <- c("ZNF683","GZMH","KLRF1","KLRD1","PDCD1","KLRC1","CX3CR1","STMN1","CCR7","LEF1","TCF7","SELL",
        "GATA3","CD40LG","CCR4",
        "ICOS","CCR5","GZMK","EOMES","LAG3","KLRB1","TBX21","CXCR3","IFNG","KLRG1","CCR6", "STAT1",
        "RORC", "IL17A","IL26","IL17F", "IL21", "FOXP3","IL2RA","CTLA4", "TNFRSF9")


sce$lvl2 <- factor(sce$lvl2, levels = c("CD4 T eff.","CD4 TEM","CD4 TCM",
                                        "CD4 TN","Th2","Th1/EMRA","Th17","Treg"))



dot_cd4 <- plotDots(sce, features = fs, group = "lvl2", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  coord_flip() + 
  theme_minimal(6) + 
  theme(
    aspect.ratio = length(unique(sce$lvl2))/length(fs),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text = element_text(size = 8, color = "black")
  ) + 
  scale_size_continuous(range = c(2, 4), breaks = c(0.25, 0.75)) 

dot_cd4

cd4 <- sce

# NK
################################################################################
################################################################################
lin <- "NK"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/{lin}/lvl2_sce.qs"), nthreads = 8)

sce$lvl2 <- as.character(sce$lvl2)
sce$lvl2[sce$lvl2 == "cycling NK"] <- "prolif. NKs"
sce$lvl2[sce$lvl2 == "CD56hi CD16- NK"] <- "CD56bright CD16- NK"
sce$lvl2[sce$lvl2 == "CD56low CD16high NK"] <- "CD56dim CD16bright NK"

sce$lvl2 <- factor(sce$lvl2, levels = c( "CD56bright CD16- NK","CD56dim CD16bright NK","prolif. NKs"))

fs <- c("NCAM1","CCR5","GZMK","KLRK1","CXCR3","KLRC1","KLRD1","KLRF1",
        "FCGR3A","KLRB1","IL12RB1","KLRG1","GZMB","GZMA","CX3CR1",
        "EOMES","CD38","STMN1","TOP2A","CD34")



dot_nk <- plotDots(sce, features = intersect(fs, rownames(sce)), group = "lvl2", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  coord_flip() + 
  theme_minimal(6) + 
  theme(
    aspect.ratio = length(unique(sce$lvl2))/length(fs),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
    panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
    axis.title = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text = element_text(size = 9, color = "black")
  ) + 
  scale_size_continuous(range = c(2, 4), breaks = c(0.25, 0.75)) 

dot_nk

nk <- sce

# DOT plots
################################################################################
################################################################################
pdf("/Users/emanuelepitino/Desktop/fig_s11/dots.pdf")
dot_cd4
dot_cd8
dot_nk
dot_b
dot_myelo
dev.off()
## UMAPS
################################################################################
################################################################################
um_func <- function(sc,lin){
cd <- as.data.frame(colData(sc))
um <- reducedDim(sc, "UMAP")
cd$um1 <- um[rownames(cd), 1]
cd$um2 <- um[rownames(cd), 2]

um <- ggplot(cd, aes(x = um1, y = um2, color = lvl2 )) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.001), dpi = 500) +
  scale_color_manual(values = pal) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "", y = "", color = glue("{lin}")) +
  theme(legend.position = "right") +
  theme(panel.grid = element_blank(),
        text = element_text(color = "black"), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01),"cm")) +
guides(color = guide_legend(ncol = 2,nrow = 5, override.aes = list(size = 2)))

return(um)
}

ctypes <- c(unique(b$lvl2),unique(myelo$lvl2),unique(cd4$lvl2),unique(cd8$lvl2), unique(as.factor(nk$lvl2)))

set.seed(12345)
pal <- Polychrome::createPalette(31,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(ctypes)

um <- wrap_plots(
  um_func(cd4,"CD4"),
  um_func(cd8,"CD8"),
  um_func(nk,"NK"),
  um_func(b,"B"),
  um_func(myelo,"Myeloid"),
  ncol = 1)

pdf("/Users/emanuelepitino/Desktop/fig_s11/umaps.pdf", width  = 5, height = 8)
um
dev.off()

















