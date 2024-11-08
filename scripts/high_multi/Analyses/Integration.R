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

sce <- qread(glue("{res_dir}/sub_25pct_merged.qs"), nthreads = 8)


sce <- harmony::RunHarmony(
  sce,
  group.by.vars = "tech",
  dims.use = 1:20,
  verbose = TRUE,
  reduction.save = "HARMONY",
  ncores = 8
)

####################################################################################
## integrated PCA plots
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
pca <- reducedDim(sce,"HARMONY") # pc

df <- cbind(df,pca) # merge
df <- df[sample(nrow(df)), ] # shuffle

pca <- \(var){
  pal2 <- Polychrome::createPalette(n_distinct(sce[[var]]), c("#FBB4AE", "#B3CDE3", "#CCEBC5"))
  names(pal2) <- unique(sce[[var]])
  
  
  p <- ggplot(df, aes(x = HARMONY_1, y = HARMONY_2, color = !!sym(var))) +
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

# Plot unique samples in PCA
unique <- names(which(tapply(sce$tech, sce$sample, function(x) length(unique(x)) == 1)))
sce$unique <- ifelse(sce$sample %in% unique,"yes","no")

df <- as.data.frame(colData(sce)) # cd
pca_df <- reducedDim(sce,"HARMONY") # pc

df <- cbind(df,pca_df) # merge
df <- df[sample(nrow(df)), ] # shuffle

df$unique <- factor(df$unique, levels = c("yes","no"))
pca_unique_samples <- ggplot(df, aes(x = HARMONY_1, y = HARMONY_2, color = unique)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 1200) +
  theme_bw() +
  scale_color_manual(values = c("yes" = "red4", "no" = "grey80")) +
  theme(
    text = element_text(size = 25, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  coord_equal() +
  labs(color = "Unique \nacross \nreplicates")




pca_plts <- wrap_plots(
  pca("replicate") + theme(legend.position = "none", axis.title = element_blank()),
  pca_unique_samples + theme(axis.title.x = element_blank()),
  pca("sample") + theme(legend.position = "none", axis.title.y = element_blank()),
  ncol = 1)

#### save
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/harmony_pca_plts.pdf"),width = 6)
pca_plts
dev.off()

 ## Lisi
#sub <- sce[,sample(colnames(sce),1000)]
unint <- reducedDim(sce,"PCA")[,1:2]
int <- reducedDim(sce,"HARMONY")[,1:2]

lisi_plt <- \(md_var){
cd <- as.data.frame(colData(sce)) %>% select(md_var)
lisi_unint <- lisi::compute_lisi(unint, cd, md_var)
lisi_int <- lisi::compute_lisi(int, cd, md_var)

df <- data.frame(var = c("PCA","HARMONY"),
                 score = c(mean(lisi_unint[[md_var]]),mean(lisi_int[[md_var]])))

df$var <- factor(df$var, levels = c("PCA","HARMONY"))
p <- ggplot(df, aes(x = var, y = score,)) +
  geom_col() +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 3/1,
        text = element_text(size = 18, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.subtitle = element_text(hjust = 1, size = 16, color = "black")) +
  labs(y = "LISI score", x = "")
if(md_var == "tech") {p <- p  + 
  scale_y_continuous(breaks = c(0.0,1,round(max(df$score),1))) +
  labs(subtitle = "Tech")
  }
if(md_var == "replicate") {p <- p  + scale_y_continuous(breaks = c(0.0,1.5,round(max(df$score),1))) +
  labs(subtitle = "Replicate")}

return(p)
}

gg_lisi <- wrap_plots(
  lisi_plt("tech"),
  lisi_plt("replicate") + theme(axis.title.y = element_blank()),
  ncol = 2) 



#### save
outdir <- glue("{plt_dir}/{stamp}")
dir.create(outdir,showWarnings = F)
pdf(glue("{outdir}/lisi_plt.pdf"),width = 6)
gg_lisi
dev.off()

qsave(sce, file = glue("{res_dir}/integrated_25pct_merged.qs"), nthreads = 8)

