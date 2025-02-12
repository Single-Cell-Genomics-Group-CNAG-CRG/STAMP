## STAMP11b + STAMP11c

stamp11b_pbmc <- readRDS('stamp11b_pbmc_seurat.rds')
stamp11c_pbmc <- readRDS('stamp11c_pbmc_seurat.rds')
stamp11b_pbmc$STAMP <- 'STAMP-XCP'
stamp11c_pbmc$STAMP <- 'STAMP-CP'

stamp1_pbmc_merged <- merge(stamp11b_pbmc, stamp11c_pbmc)
stamp1_pbmc_merged <- JoinLayers(stamp1_pbmc_merged)

## background merged
bkg_stamp11_merged <- rbind(bkg_stamp11b, bkg_stamp11c)

## Label low quality cells

stamp1_pbmc_merged$Quality <- 'Good Quality'
stamp1_pbmc_merged$Quality[stamp1_pbmc_merged$nCount_CP > 7000 | stamp1_pbmc_merged$nCount_CP < 200 |
                             stamp1_pbmc_merged$Area.um2 < 40 | stamp1_pbmc_merged$Area.um2 > 250] <- 'Low Quality'
table(stamp1_pbmc_merged$Quality, stamp1_pbmc_merged$STAMP)



#QC Plots
#
# VlnPlot(stamp1_pbmc_merged, features = c('nCount_CP', 'Area.um2', 'Mean.PanCK', 'Mean.DAPI'), group.by = 'STAMP', pt.size = 0, ncol = 2) & scale_y_log10() &
#   scale_fill_manual(values = c('#898999', '#C6431F')) &  theme(text=element_text(size=16))

#pdf('/mnt/scratch2/STAMP/Figures/STAMP11_QC_Violins.pdf', 8,4)
p <- VlnPlot(stamp1_pbmc_merged, features = c('nCount_CP', 'Area.um2'), group.by = 'STAMP', pt.size = 0, ncol = 2) & scale_y_log10() &
  scale_fill_manual(values = c('#898999', '#C6431F')) &  theme(text=element_text(size=16))
#dev.off()
ggsave(file="/mnt/scratch2/STAMP/Figures/STAMP11_QC_Violins.svg", plot=p, width=8, height=4)




## Barplot of cell numbers
pdf('/mnt/scratch2/STAMP/Figures/STAMP11_QC_Barplot_NCells.pdf', 5,8)
ggplot(stamp1_pbmc_merged@meta.data, aes(STAMP, fill = Quality)) +
  geom_bar(position = "fill") +            # Use position = "fill" for proportions
  ggtitle("Cell quality - CP") +         # Update title to reflect proportions
  scale_fill_manual(values = c('#89B790', '#05485E')) +
  theme_bw() + ylab('Proportion') +
  theme(text = element_text(size = 16))
dev.off()



## Plot mean expression of every protein per Run (STAMp11b , STAMP11c)
average_expression <- AverageExpression(stamp1_pbmc_merged, group.by = "STAMP")$CP

# Convert to a data frame
avg_expr_df <- as.data.frame((average_expression))


#pdf('/mnt/scratch2/STAMP/Figures/STAMP11_CorrelationAvgExpression.pdf', 8,4)
p1 <- ggplot(avg_expr_df, aes(x = `STAMP-CP`, y = `STAMP-XCP`)) +
  geom_point(color = 'gray', size = 5) +
  geom_smooth(method = "lm", color = "red", se = F, linetype = 'dashed') +
  stat_cor(method = "pearson", label.x = 3, label.y = 1000) +
  labs(x = "Average Expression (STAMP-CP)", y = "Average Expression (STAMP-X-CP)", title = "") +
   theme_bw() +  theme(text=element_text(size=16)) +
  theme(legend.position = "none")  + coord_fixed() + scale_x_continuous( limits = c(0,2000)) + scale_y_continuous(limits = c(0, 2000))#Remove legend if coloring by gene names
#dev.off()
ggsave(file="/mnt/scratch2/STAMP/Figures/STAMP11_CorrelationAvgExpression.png", plot=p, width=6, height=4)



stamp1_pbmc_merged <- NormalizeData(stamp1_pbmc_merged)
stamp1_pbmc_merged <- ScaleData(stamp1_pbmc_merged)
stamp1_pbmc_merged <- FindVariableFeatures(stamp1_pbmc_merged)

stamp1_pbmc_merged <- SketchData(
  object = stamp1_pbmc_merged,
  ncells = 100000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(stamp1_pbmc_merged) <- "sketch"
stamp1_pbmc_merged <- FindVariableFeatures(stamp1_pbmc_merged)
stamp1_pbmc_merged <- ScaleData(stamp1_pbmc_merged)
stamp1_pbmc_merged <- RunPCA(stamp1_pbmc_merged)
stamp1_pbmc_merged <- FindNeighbors(stamp1_pbmc_merged, dims = 1:15)
stamp1_pbmc_merged <- FindClusters(stamp1_pbmc_merged, resolution = 0.5)
stamp1_pbmc_merged <- RunUMAP(stamp1_pbmc_merged, dims = 1:15, return.model = T)

library(harmony)
stamp1_pbmc_merged <- harmony::RunHarmony(stamp1_pbmc_merged, 'STAMP')
stamp1_pbmc_merged <- FindNeighbors(stamp1_pbmc_merged, dims = 1:15, reduction = 'harmony')
stamp1_pbmc_merged <- FindClusters(stamp1_pbmc_merged, resolution = 0.5)
stamp1_pbmc_merged <- RunUMAP(stamp1_pbmc_merged, dims = 1:15, return.model = T, reduction = 'harmony')

stamp1_pbmc_merged <- ProjectData(
  object = stamp1_pbmc_merged,
  assay = "CP",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)
DefaultAssay(stamp1_pbmc_merged) <- 'CP'


DimPlot(stamp1_pbmc_merged, label = F, group.by = c('STAMP', 'seurat_clusters'))
DimPlot(stamp1_pbmc_merged, label = F, group.by = c('STAMP', 'cluster_full'), reduction = 'full.umap')


dittoBarPlot(stamp1_pbmc_merged, var = 'cluster_full', group.by = 'STAMP')


##### CEll type annotation
# match STAMPX (RNA) - nomenclature "B"       "CD4 T"   "CD8 T"   "DC"      "Mono"    "NK"


DotPlot(stamp1_pbmc_merged, features = rownames(stamp1_pbmc_merged), scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

stamp1_pbmc_merged$CellType_pred <- as.character(stamp1_pbmc_merged$cluster_full)

stamp1_pbmc_merged[[]] <- stamp1_pbmc_merged[[]] %>%
  mutate(CellType_pred = case_when(CellType_pred %in% c('0') ~ 'CD4 T',
                                   CellType_pred %in% c('3','4','10') ~ 'CD8 T',
                                   CellType_pred %in% '2' ~ 'B',
                                   CellType_pred %in% c('1', '9') ~ 'Mono',
                                   CellType_pred %in% '5' ~ 'NK',
                                   CellType_pred %in% c('6','7', '8') ~ 'other T',
                                   CellType_pred %in% '11' ~ 'DC',
                                   CellType_pred %in% '12' ~ 'Unresolved'))

saveRDS(stamp1_pbmc_merged, file = 'STAMP11_b_c_merged_annotated.rds')


DotPlot(stamp1_pbmc_merged, features = rownames(stamp1_pbmc_merged), group.by = 'CellType_pred', scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()




 ## UMAP plots for figure
p1 <- DimPlot(stamp1_pbmc_merged, group.by = c('cluster_full'), label = F,reduction = 'full.umap', label.size = 5) + ggtitle('') & NoAxes()
p2 <- DimPlot(stamp1_pbmc_merged, group.by = c('STAMP'), label = F,reduction = 'full.umap', cols = c('#898999', '#C6431F')) + ggtitle('') & NoAxes()
p3 <- DimPlot(stamp1_pbmc_merged, group.by = c('CellType_pred'), label = F,reduction = 'full.umap', cols = brewer.pal(8, 'Paired')) + ggtitle('') & NoAxes()

png('/mnt/scratch2/STAMP/Figures/STAMP11_UMAPs.png', width = 8, height = 4, units = 'in', res = 300)
p2  + p3
dev.off()


## Barplot cell proportions

proportions_df <- stamp1_pbmc_merged[[]] %>%
  group_by(STAMP, cluster_full) %>%
  summarise(count = n()) %>%
  group_by(STAMP) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Step 2: Create the barplot
p4 <- ggplot(proportions_df, aes(x = STAMP, y = proportion, fill = cluster_full)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "STAMP", y = "Proportion", fill = "Cluster") +
  theme_bw() + theme(text=element_text(size=16)) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   & NoLegend()

proportions_df_cellt <- stamp1_pbmc_merged[[]] %>%
  group_by(STAMP, CellType_pred) %>%
  summarise(count = n()) %>%
  group_by(STAMP) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p5 <- ggplot(proportions_df_cellt, aes(x = STAMP, y = proportion, fill = CellType_pred)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "STAMP", y = "Proportion", fill = "CellType") +
  theme_bw() +theme(text=element_text(size=16)) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + scale_fill_manual(values = brewer.pal(8, 'Paired')) & NoLegend()

library(patchwork)


Idents(stamp1_pbmc_merged) <- stamp1_pbmc_merged$CellType_pred
markers <- FindAllMarkers(stamp1_pbmc_merged)
markers_11b <- FindAllMarkers(subset(stamp1_pbmc_merged, STAMP == 'STAMP11b'))
markers_11c <- FindAllMarkers(subset(stamp1_pbmc_merged, STAMP == 'STAMP11c'))


markers_11b <- markers_11b[markers_11b$avg_log2FC>0, ]
markers_11c <- markers_11c[markers_11c$avg_log2FC>0, ]

top_markers <- markers %>% group_by(cluster) %>% top_n(4, avg_log2FC) %>% as.data.frame()
top_markers_11b <- markers_11b %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% as.data.frame()
top_markers_11c <- markers_11c %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% as.data.frame()

DotPlot(stamp1_pbmc_merged, features = unique(top_markers$gene), group.by = 'CellType_pred', scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()






# DotPlot figure
DotPlot(stamp1_pbmc_merged, features = c('CD19','CD20','CD4', 'CD3', 'CD27', 'CD45',
                                         'CD8' ,'CCR7', 'Vimentin', 'TCF7',
                                         'CD11b', 'CD11c', 'CD14', 'CD16', 'CD38',
                                         'CD56', 'CD45RA', 'pan-RAS'), group.by = 'CellType_pred') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis() + coord_flip()




mat_hm <- stamp1_pbmc_merged@assays$CP$data[c('CD19','CD20','CD4', 'CD3', 'CD27', 'CD45',
                                              'CD8' ,'CCR7', 'Vimentin', 'TCF7',
                                              'CD11b', 'CD11c', 'CD14', 'CD16', 'CD38',
                                              'CD56', 'CD45RA', 'pan-RAS'), ]


pheatmap(mat_hm, border_color = NA, cluster_rows = F, cluster_cols = F, show_colnames = F)





avg_stamp11 <- AverageExpression(stamp1_pbmc_merged, group.by = "CellType_pred", layer = 'data')$CP
avg_stamp11_clust <- AverageExpression(stamp1_pbmc_merged, group.by = "cluster_full", layer = 'data')$CP


color = colorRampPalette(c("navy", "white", "red"))(21)
breaks = seq(-2,2, 0.2)


png('/mnt/scratch2/STAMP/Figures/STAMP11_CellType_HM.png', width = 8, height = 8, units = 'in', res = 300)
pheatmap(avg_stamp11[c('CD19','CD20','CD4', 'CD3', 'CD27', 'CD45',
                       'CD8' ,'CCR7', 'Vimentin',
                       'CD11c', 'CD14', 'CD16', 'CD38',
                       'CD45RA', 'pan-RAS'),], scale = 'row', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = F, fontsize = 16, angle_col = 45)
dev.off()




pheatmap(avg_stamp11, scale = 'row', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = F, fontsize = 16, angle_col = 45)
pheatmap(avg_stamp11_clust, scale = 'row', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = F, fontsize = 16, angle_col = 45)









