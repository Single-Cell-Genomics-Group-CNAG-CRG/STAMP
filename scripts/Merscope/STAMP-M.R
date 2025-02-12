### STAMP M
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)

path = '/mnt/scratch2/STAMP/STAMP-M/'
## There are 9 folders/output for 12 stamps.
samples <- dir(path)

stampm <- list()
for(i in samples) {
  stampm[[i]]$exp <- fread(paste0(path, i, '/cell_by_gene.csv'))
  # get cell id
  cell_id <- as.character(stampm[[i]]$exp$cell)
  # transpose
  stampm[[i]]$exp <- t( stampm[[i]]$exp[, -1])
  colnames(stampm[[i]]$exp) <- cell_id
  ## remove blank probes
  w <- (grep('Blank', rownames(stampm[[i]]$exp)))
  stampm[[i]]$exp <- stampm[[i]]$exp[-w, ]

  stampm[[i]]$meta <- as.data.frame(fread(paste0(path, i, '/cell_metadata.csv')))
  rownames( stampm[[i]]$meta) <- as.character( stampm[[i]]$meta$EntityID)

  stampm[[i]]<- CreateSeuratObject(counts = stampm[[i]]$exp, assay = 'RNA',
                                                 meta.data =stampm[[i]]$meta )
  # add name
  stampm[[i]]$sample <- i

}

# region_0 needs to be separated into 3 objects (3 stamps)
stampm$region_0$STAMPs <- "region_0_left"
stampm$region_0$STAMPs[stampm$region_0$center_x > 6400 & stampm$region_0$center_x < 9850] <- 'region_0_mid'
stampm$region_0$STAMPs[stampm$region_0$center_x > 9850 ] <- 'region_0_right'
region0 <- SplitObject(stampm$region_0, split.by = 'STAMPs')

## remove original region 0 and add new list to existing list

stampm <- stampm[-1]
stampm <- c(region0, stampm)

## rename based on assumption
#names(stampm) <- c('L', 'K', 'Q', 'J', 'M', 'I', 'P', 'S', 'MX1', 'D', 'C', 'B')
# merge into one
all <- Reduce(function(x, y) merge(x, y), stampm)
all <- JoinLayers(all)
all$sample[all$sample == 'region_0'] <- all$STAMPs[all$sample == 'region_0']

all$cellline <- as.factor(all$sample)
#levels(all$sample) <- c('L', 'K', 'Q', 'J', 'M', 'I', 'P', 'S', 'MX1', 'D', 'C', 'B')
levels(all$cellline) <- c('L', 'K', 'Q', 'J', 'S', 'P', 'I', 'M', 'MX1', 'D', 'C', 'B')
all$cellline <- as.character(all$cellline)
# area
all$area <- (all$max_x - all$min_x) * (all$max_y - all$min_y)


p <-VlnPlot(all, features = c('nCount_RNA', 'nFeature_RNA', 'area'), pt.size = 0, group.by = 'cellline') &
  stat_summary(fun = median, geom = "point", size = 2, color = "black") &
  scale_fill_manual(values = as.character(Polychrome::kelly.colors(13)[2:13]))
p <- p & scale_y_log10()

# add qc flag
all$qcFlag <- 'High quality'
all$qcFlag[all$nCount_RNA < 20 | all$nFeature_RNA < 20 | all$area > 800 | all$area < 5 ] <- 'Low quality'


summary_stats <- all[[]] %>%
  group_by(cellline) %>%
  summarize(
    n_cells = n(),  # Number of cells per group
    min_nCount_RNA = min(nCount_RNA),
    max_nCount_RNA = max(nCount_RNA),
    median_nCount_RNA = median(nCount_RNA),
    min_nFeature_RNA = min(nFeature_RNA),
    max_nFeature_RNA = max(nFeature_RNA),
    median_nFeature_RNA = median(nFeature_RNA)
  )

hist(all$nCount_RNA[all$qcFlag == 'High quality'], breaks = 50)
hist(all$nFeature_RNA, breaks = 50)


all_sub <- subset(all, qcFlag == 'High quality')
all_sub <- NormalizeData(all_sub)

Idents(all_sub) <- all_sub$cellline
marker_lines <- FindAllMarkers(all_sub)
top_markers <- marker_lines %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) %>% ungroup()

p4 <- DotPlot(all_sub, features = unique(top_markers$gene)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

p1 <- ggplot(all_sub@meta.data, aes(center_x, -center_y, color = cellline)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = as.character(Polychrome::kelly.colors(13)[2:13])) +
  theme_bw() +
  coord_fixed() +
  theme(
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.title = element_text(size = 14),  # Adjust legend title size
    axis.text.x = element_blank(),           # Remove x-axis tick labels
    axis.text.y = element_blank(),           # Remove y-axis tick labels
    axis.ticks.x = element_blank(),          # Remove x-axis ticks
    axis.ticks.y = element_blank()           # Remove y-axis ticks
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(size = 4)   # Increase point size in the legend
    )
  ) + xlab('x') + ylab('y')




p2 <-  ggplot(all_sub@meta.data , aes(center_x, -center_y, col = log10(nCount_RNA + .1))) + geom_point(size = .2)+
  theme_bw() + scale_color_viridis_c() + coord_fixed() +
  theme(
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.title = element_text(size = 14),  # Adjust legend title size
    axis.text.x = element_blank(),           # Remove x-axis tick labels
    axis.text.y = element_blank(),           # Remove y-axis tick labels
    axis.ticks.x = element_blank(),          # Remove x-axis ticks
    axis.ticks.y = element_blank()           # Remove y-axis ticks
  ) + xlab('x') + ylab('y')


p3 <-  ggplot(all_sub@meta.data , aes(center_x, -center_y, col = nFeature_RNA)) + geom_point(size = .2)+
  theme_bw() + scale_color_viridis_c(option = 'B')+ coord_fixed() +
  theme(
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.title = element_text(size = 14),  # Adjust legend title size
    axis.text.x = element_blank(),           # Remove x-axis tick labels
    axis.text.y = element_blank(),           # Remove y-axis tick labels
    axis.ticks.x = element_blank(),          # Remove x-axis ticks
    axis.ticks.y = element_blank()           # Remove y-axis ticks
  ) + xlab('x') + ylab('y')

library(patchwork)
combined_plot <- (p1 / p2 / p3) | (p / p4)

# Customize the layout for A4 format
combined_plot <- combined_plot +
  plot_layout(
    widths = c(0.6, 1.4),  # Equal width for both columns
    heights = c(1, 1, 1)  # Equal height for rows in the first column
  )   + plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = "")&
  theme(plot.tag = element_text(face = "plain")) &theme(plot.tag = element_text(size = 20))

ggsave(plot = combined_plot, filename  = 'STAMP-M_figure.png', width = 20, height = 10, dpi = 300 )





