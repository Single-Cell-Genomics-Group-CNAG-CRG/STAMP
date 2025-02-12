### STAMP 12c
source('SNR.R')
stamp11_pcf <- data.table::fread('/mnt/scratch2/STAMP/STAMP_PCF/24066_segmentation.tsv')

ggplot(stamp11_pcf, aes(`Centroid X µm`,`Centroid Y µm` )) + geom_point(size = 0.01) + geom_hline(yintercept = 1000, linetype = 'dashed') + geom_hline(yintercept = 12500, linetype = 'dashed') +
  geom_vline(xintercept = 1250, linetype = 'dashed') + geom_vline(xintercept = 9000,  linetype = 'dashed')


stamp11_pcf_pbmc <- stamp11_pcf[stamp11_pcf$`Centroid X µm` > 1250 & stamp11_pcf$`Centroid X µm` < 9000 & stamp11_pcf$`Centroid Y µm`> 1000 & stamp11_pcf$`Centroid Y µm` <12500 ,]
stamp11_pcf_pbmc <- as.data.frame(stamp11_pcf_pbmc)

ggplot(stamp11_pcf_pbmc, aes(`Centroid X µm`,`Centroid Y µm` )) + geom_point(size = 0.01) + geom_hline(yintercept = 1000, linetype = 'dashed') + geom_hline(yintercept = 12500, linetype = 'dashed') +
  geom_vline(xintercept = 1250, linetype = 'dashed') + geom_vline(xintercept = 9000,  linetype = 'dashed')


w <- grep('Cell: Mean', colnames(stamp11_pcf_pbmc))
stamp11_pcf_pbmc <- stamp11_pcf_pbmc[, w]
colnames(stamp11_pcf_pbmc) <- gsub(': Cell: Mean', '', colnames(stamp11_pcf_pbmc))
# remove DAPI
stamp11_pcf_pbmc <- stamp11_pcf_pbmc[, -1]

# add cell id
rownames(stamp11_pcf_pbmc) <- paste0('Cell', rownames(stamp11_pcf_pbmc))


stamp11_pcf_pbmc <- CreateSeuratObject(counts = Matrix::as.matrix(t(stamp11_pcf_pbmc), sparse = T), assay = 'PCF')
#Calculate SNR
snr_stamp12c <- getSNR(stamp11_pcf_pbmc)
snr_stamp12c$STAMP <- 'STAMP-PhF'

# metadata
meta <- stamp11_pcf[stamp11_pcf$`Centroid X µm` > 1250 & stamp11_pcf$`Centroid X µm` < 9000 & stamp11_pcf$`Centroid Y µm`> 1000 & stamp11_pcf$`Centroid Y µm` <12500 , 8:17]
stamp11_pcf_pbmc$Area<- meta$`Cell: Area µm^2`

VlnPlot(stamp11_pcf_pbmc, features = c('Area.um2', 'nCount_PCF', 'nFeature_PCF'), pt.size = 0)

stamp11_pcf_pbmc <- subset(stamp11_pcf_pbmc, Area.um2 > 20 &  Area.um2 < 250  & nCount_PCF > 200 & nCount_PCF < 2000)

stamp11_pcf_pbmc <- NormalizeData(stamp11_pcf_pbmc, normalization.method = 'CLR', margin = 1)
stamp11_pcf_pbmc <- ScaleData(stamp11_pcf_pbmc, do.center = T)
stamp11_pcf_pbmc <- FindVariableFeatures(stamp11_pcf_pbmc)

stamp11_pcf_pbmc <- SketchData(
  object = stamp11_pcf_pbmc,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(stamp11_pcf_pbmc) <- "sketch"
stamp11_pcf_pbmc <- FindVariableFeatures(stamp11_pcf_pbmc)
stamp11_pcf_pbmc <- ScaleData(stamp11_pcf_pbmc, do.center = T)
stamp11_pcf_pbmc <- RunPCA(stamp11_pcf_pbmc)
stamp11_pcf_pbmc <- FindNeighbors(stamp11_pcf_pbmc, dims = 1:15)
stamp11_pcf_pbmc <- FindClusters(stamp11_pcf_pbmc, resolution = 0.5)
stamp11_pcf_pbmc <- RunUMAP(stamp11_pcf_pbmc, dims = 1:15, return.model = T)

DimPlot(stamp11_pcf_pbmc)

stamp11_pcf_pbmc <- ProjectData(
  object = stamp11_pcf_pbmc,
  assay = "PCF",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)
DefaultAssay(stamp11_pcf_pbmc) <- 'PCF'

DimPlot(stamp11_pcf_pbmc, label = T, group.by = 'cluster_full', raster = F) + ggtitle( 'STAMP-PhF')

DotPlot(stamp11_pcf_pbmc, features = rownames(stamp11_pcf_pbmc), ) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

# Cell typ[e manual
stamp11_pcf_pbmc$CellType_pred <- stamp11_pcf_pbmc$cluster_full
stamp11_pcf_pbmc[[]] <- stamp11_pcf_pbmc[[]] %>%
  mutate(CellType_pred = case_when(CellType_pred %in% c('0', '4') ~ 'CD8 T',
                                   CellType_pred %in% c('1','7', '17') ~ 'CD4 T',
                                   CellType_pred %in% c('2', '15', '3', '13') ~ 'Mono',
                                   CellType_pred %in% '5' ~ 'B',
                                   CellType_pred %in% c('9', '10', '14') ~ 'NK',
                                   CellType_pred %in% c('6', '8', '11') ~ 'other T',
                                   CellType_pred %in% '11' ~ 'DC',
                                   CellType_pred %in% c('12','16', '18', '19') ~ 'other'))

DefaultAssay(stamp11_pcf_pbmc) <- 'PCF'
DimPlot(stamp11_pcf_pbmc, label = T, group.by = c('cluster_full', 'CellType_pred'))


#pdf('/mnt/scratch2/STAMP/Figures/STAMP1_PCF_CTCIdentification.pdf', 8,6)
png('/mnt/scratch2/STAMP/Figures/STAMP1_PCF_CTCIdentification..png', width = 8, height = 6, units = 'in', res = 300)
FeatureScatter(stamp11_pcf_pbmc, feature1 = 'Area.um2', feature2 = 'Pan-Cytokeratin', group.by = 'orig.ident', cols = 'grey', raster = F) +
  ggtitle('CTC identification') + geom_hline(yintercept = 7, linetype = 'dashed', color = 'red') + geom_vline(xintercept = 210, linetype = 'dashed', color = 'red')  +
  theme_bw() +theme(text=element_text(size=16)) & NoLegend()
dev.off()



## STAMP 12b
stamp12b <- data.table::fread('/mnt/scratch2/STAMP/STAMP_PCF_Qupath/stamp12b_segmentation.tsv')


ggplot(stamp12b, aes(`Centroid X µm`,`Centroid Y µm` )) + geom_point(size = 0.01) +
  geom_hline(yintercept = 1000, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 12000, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 1250, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 8000,  linetype = 'dashed', color = 'red')


stamp12b_pcf <- stamp12b[stamp12b$`Centroid X µm` > 1250 & stamp12b$`Centroid X µm` < 8000 & stamp12b$`Centroid Y µm`> 1000 & stamp12b$`Centroid Y µm` <12000 ,]
stamp12b_pcf <- as.data.frame(stamp12b_pcf)

meta_stamp12b <- stamp12b[stamp12b$`Centroid X µm` > 1250 & stamp12b$`Centroid X µm` < 8000 & stamp12b$`Centroid Y µm`> 1000 & stamp12b$`Centroid Y µm` <12000 , 8:17]


w <- grep('Cell: Mean', colnames(stamp12b_pcf))
stamp12b_pcf <- stamp12b_pcf[, w]
colnames(stamp12b_pcf) <- gsub(': Cell: Mean', '', colnames(stamp12b_pcf))
# remove DAPI
stamp12b_pcf <- stamp12b_pcf[, -1]

# add cell id
rownames(stamp12b_pcf) <- paste0('Cell', rownames(stamp12b_pcf))


stamp12b_pcf_pbmc <- CreateSeuratObject(counts = Matrix::as.matrix(t(stamp12b_pcf), sparse = T), assay = 'PCF')
stamp12b_pcf_pbmc$Area <- meta_stamp12b$`Cell: Area µm^2`
# Calculate SNR
snr_stamp12b <- getSNR(stamp12b_pcf_pbmc)
snr_stamp12b$STAMP <- 'STAMP-XPhF'

VlnPlot(stamp12b_pcf_pbmc, features = c('nCount_STAMP12b', "Area" ), pt.size = 0)

# Filter
stamp12b_pcf_pbmc <- subset(stamp12b_pcf_pbmc, Area > 20 &  Area < 250  & nCount_PCF > 200 & nCount_PCF < 2000)

# Workflow
stamp12b_pcf_pbmc <- NormalizeData(stamp12b_pcf_pbmc)
stamp12b_pcf_pbmc <- ScaleData(stamp12b_pcf_pbmc)
stamp12b_pcf_pbmc <- FindVariableFeatures(stamp12b_pcf_pbmc)

stamp12b_pcf_pbmc <- SketchData(
  object = stamp12b_pcf_pbmc,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(stamp12b_pcf_pbmc) <- "sketch"
stamp12b_pcf_pbmc <- FindVariableFeatures(stamp12b_pcf_pbmc)
stamp12b_pcf_pbmc <- ScaleData(stamp12b_pcf_pbmc)
stamp12b_pcf_pbmc <- RunPCA(stamp12b_pcf_pbmc)
stamp12b_pcf_pbmc <- FindNeighbors(stamp12b_pcf_pbmc, dims = 1:15)
stamp12b_pcf_pbmc <- FindClusters(stamp12b_pcf_pbmc, resolution = 0.5)
stamp12b_pcf_pbmc <- RunUMAP(stamp12b_pcf_pbmc, dims = 1:15, return.model = T)

DimPlot(stamp12b_pcf_pbmc)

stamp12b_pcf_pbmc <- ProjectData(
  object = stamp12b_pcf_pbmc,
  assay = "PCF",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "seurat_clusters")
)
DefaultAssay(stamp12b_pcf_pbmc) <- 'PCF'

DimPlot(stamp12b_pcf_pbmc, label = T, group.by = 'cluster_full', raster = F) + ggtitle('STAMP-XPhF')

DotPlot(stamp12b_pcf_pbmc, features = unique(rownames(stamp12b_pcf_pbmc)) ) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()



## Average expression per cluster


avg_stamp12b <- AverageExpression(stamp12b_pcf_pbmc, group.by = "cluster_full", layer = 'data')$STAMP12b
avg_stamp12c <- AverageExpression(stamp11_pcf_pbmc, group.by = "cluster_full", layer = 'data')$PCF

p1 <- pheatmap(avg_stamp12b, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)

p2 <- pheatmap(avg_stamp12c, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)

gridExtra::grid.arrange(p1$gtable ,p2$gtable, ncol = 2)



## Merged

stamp12b_pcf_pbmc$Sample = 'STAMP-XPhF'
stamp11_pcf_pbmc$Sample = 'STAMP-PhF'

stamp12_merged <- merge(stamp12b_pcf_pbmc, stamp11_pcf_pbmc)
stamp12_merged <- JoinLayers(stamp12_merged)
stamp12_merged <- NormalizeData(stamp12_merged, normalization.method = 'CLR', margin = 1)
#rownames(stamp12_merged@meta.data) <- Cells(stamp12_merged)


stamp12_merged$Area.um2[stamp12_merged$Sample %in% 'STAMP-XPhF'] <- stamp12_merged$Area[stamp12_merged$Sample %in% 'STAMP-XPhF']

#pdf('/mnt/scratch2/STAMP/Figures/STAMP12_QC_Violins.pdf', 8,4)
p <- VlnPlot(stamp12_merged, features = c('nCount_PCF', 'Area.um2'), group.by = 'Sample', pt.size = 0, ncol = 2) & scale_y_log10() &
  scale_fill_manual(values = c('#898999', '#C6431F')) &  theme(text=element_text(size=16))
#dev.off()
ggsave(file="/mnt/scratch2/STAMP/Figures/STAMP12_QC_Violins.svg", plot=p, width=8, height=4)



stamp12_merged <- stamp12_merged %>% NormalizeData() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:15)

DimPlot(stamp12_merged, group.by = 'Sample')
DimPlot(stamp12_merged, group.by = 'seurat_clusters')

library(harmony)
stamp12_merged <- harmony::RunHarmony(stamp12_merged, 'Sample')
stamp12_merged <- FindNeighbors(stamp12_merged, dims = 1:15, reduction = 'harmony')
stamp12_merged <- FindClusters(stamp12_merged, resolution = 0.5)
stamp12_merged <- RunUMAP(stamp12_merged, dims = 1:15, return.model = T, reduction = 'harmony')
DimPlot(stamp12_merged, group.by = c('Sample', 'seurat_clusters'))

write_rds(stamp12_merged, file = 'stamp12_pcf_merged.rds')


dittoSeq::dittoBarPlot(stamp12_merged, var = 'seurat_clusters', group.by = 'Sample')

## Plot mean expression of every protein per Run (STAMp11b , STAMP11c)
average_expression <- AverageExpression(stamp12_merged, group.by = "Sample")$PCF

# Convert to a data frame
avg_expr_df <- as.data.frame((average_expression))

#pdf('/mnt/scratch2/STAMP/Figures/STAMP12_CorrelationAvgExpression.pdf', 8,4)
p2 <- ggplot(avg_expr_df, aes(x = `STAMP-PhF`, y = `STAMP-XPhF`)) +
  geom_point(color = 'gray', size = 5) +
  geom_smooth(method = "lm", color = "red", se = F, linetype = 'dashed') +
  stat_cor(method = "pearson", label.x = 3, label.y = 1000) +
  labs(x = "Average Expression (STAMP-PhF)", y = "Average Expression (STAMP-XPhF)",
       title = "") +
  theme_bw() +  theme(text=element_text(size=16)) +
  theme(legend.position = "none")  + coord_fixed()  + scale_x_continuous( limits = c(0,2000)) + scale_y_continuous(limits = c(0, 2000))# Remove legend if coloring by gene names
#dev.off()
#ggsave(file="/mnt/scratch2/STAMP/Figures/STAMP12_CorrelationAvgExpression.png", plot=p, width=6, height=4)
ggsave(file="/mnt/scratch2/STAMP/Figures/STAMP_CorrelationAvgExpression_CP_PhF.png", plot=p1+p2, width=10, height=6)


avg_stamp12_merged <- AverageExpression(stamp12_merged, group.by = "cluster_full", layer = 'data')$PCF
pheatmap(avg_stamp12_merged, scale = 'row', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)


## Barplot proportion of low quality cells per sample
# have to remake the merged object to label low quality cells

stamp12b_raw <- stamp12b[stamp12b$`Centroid X µm` > 1250 & stamp12b$`Centroid X µm` < 8000 & stamp12b$`Centroid Y µm`> 1000 & stamp12b$`Centroid Y µm` <12000 ,]
stamp12b_raw <- as.data.frame(stamp12b_raw)
meta_stamp12b <- stamp12b[stamp12b$`Centroid X µm` > 1250 & stamp12b$`Centroid X µm` < 8000 & stamp12b$`Centroid Y µm`> 1000 & stamp12b$`Centroid Y µm` <12000 , 8:17]
w <- grep('Cell: Mean', colnames(stamp12b_raw))
stamp12b_raw <- stamp12b_raw[, w]
colnames(stamp12b_raw) <- gsub(': Cell: Mean', '', colnames(stamp12b_raw))
# remove DAPI
stamp12b_raw <- stamp12b_raw[, -1]
# add cell id
rownames(stamp12b_raw) <- paste0('Cell', rownames(stamp12b_raw))
stamp12b_raw <- CreateSeuratObject(counts = Matrix::as.matrix(t(stamp12b_raw), sparse = T), assay = 'PCF')
stamp12b_raw$Area <- meta_stamp12b$`Cell: Area µm^2`
stamp12b_raw$STAMP = 'STAMP-XPhF'

# STAMP12c
stamp12c_raw <- stamp11_pcf[stamp11_pcf$`Centroid X µm` > 1250 & stamp11_pcf$`Centroid X µm` < 9000 & stamp11_pcf$`Centroid Y µm`> 1000 & stamp11_pcf$`Centroid Y µm` <12500 ,]
stamp12c_raw <- as.data.frame(stamp12c_raw)
meta <- stamp11_pcf[stamp11_pcf$`Centroid X µm` > 1250 & stamp11_pcf$`Centroid X µm` < 9000 & stamp11_pcf$`Centroid Y µm`> 1000 & stamp11_pcf$`Centroid Y µm` <12500 , 8:17]
w <- grep('Cell: Mean', colnames(stamp12c_raw))
stamp12c_raw <- stamp12c_raw[, w]
colnames(stamp12c_raw) <- gsub(': Cell: Mean', '', colnames(stamp12c_raw))
# remove DAPI
stamp12c_raw <- stamp12c_raw[, -1]
# add cell id
rownames(stamp12c_raw) <- paste0('Cell', rownames(stamp12c_raw))
stamp12c_raw <- CreateSeuratObject(counts = Matrix::as.matrix(t(stamp12c_raw), sparse = T), assay = 'PCF')
stamp12c_raw$Area <- meta$`Cell: Area µm^2`
stamp12c_raw$STAMP <- 'STAMP-PhF'

avg_stamp12c <- AverageExpression(stamp11_pcf_pbmc, group.by = 'orig.ident')
write.csv(avg_stamp12c, quote = F, file = 'avg_stamp12c.csv')

stamp12_raw_merged <- merge(stamp12b_raw, stamp12c_raw)

stamp12_raw_merged$Quality <- 'Good Quality'
stamp12_raw_merged$Quality[stamp12_raw_merged$nCount_PCF > 2000 |  stamp12_raw_merged$nCount_PCF < 200 | stamp12_raw_merged$Area < 20 |
                             stamp12_raw_merged$Area > 250] <- 'Low Quality'

ggplot(stamp12_raw_merged@meta.data, aes(STAMP, fill = Quality)) +
  geom_bar(position = "fill") +            # Use position = "fill" for proportions
  ggtitle("Cell quality - STAMP PhF ") +         # Update title to reflect proportions
  scale_fill_manual(values = c('#89B790', '#05485E')) +
  theme_bw() + ylab('Proportion') +
  theme(text = element_text(size = 16))




### SNR for both samples - visualization ##
snr_all <- data.frame(row.names = rownames(snr_stamp12b), `SNR_STAMP-XPhf` = snr_stamp12c$Top_Bottom_Ratio,
                      `SNR_STAMP-Phf` = snr_stamp12b$Top_Bottom_Ratio )

snr_all <- pivot_longer(snr_all, cols = everything(), names_to = "STAMP", values_to = "SNR")



# Create the boxplot with connected lines
ggplot(snr_long, aes(x = STAMP, y = Value, group = Condition, color = Condition)) +
  geom_boxplot() +
  geom_line(aes(group = Condition), position = position_dodge(width = 0.75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Gene", y = "Value", title = "Boxplot with Connected Lines Between Points")


















