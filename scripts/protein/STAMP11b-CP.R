### STAMP 11b
library(tidyverse)
library(Seurat)
library(viridis)
library(plotly)
library(data.table)

## Read in seurat manually

cp_counts <- data.table::fread('/mnt/scratch2/STAMP/stamp_11b/STAMPXCP_exprMat_file.csv.gz')
cp_counts <- cp_counts[, -c(1:2)] # remove fov and cell id
aux_cp_counts <- cp_counts[, c(31:35,68:69)]
cp_counts <- cp_counts[, -c(31:35,68:69)] # remove membrane markers and control iGGs
cp_meta <- data.table::fread('/mnt/scratch2/STAMP/stamp_11b/STAMPXCP_metadata_file.csv.gz')
cp_counts <- as.data.frame(t(cp_counts))
colnames(cp_counts) <- cp_meta$cell_id
cell_meta <- data.table::fread('/mnt/scratch2/STAMP/stamp_11b/STAMPXCP-polygons.csv.gz') %>% as.data.frame
cell_meta$cell <- make.names(cell_meta$cell, unique = T)
cell_meta <- cell_meta[cell_meta$cell %in% colnames(cp_counts) , ]
aux_cp_counts <- as.data.frame(t(aux_cp_counts))
colnames(aux_cp_counts) <- cp_meta$cell_id

#Make data frame of only background
bkg_stamp11b <- data.frame('MsIgG1_STAMP11b' = t(aux_cp_counts[6,]),
                           'RbIgG_STAMP11b' = t(aux_cp_counts[7, ]))

bkg_stamp11b$STAMP <- 'STAMP11b'

## Remove background from counts (last 2 rows of aux_cp_counts)
cp_counts <- cp_counts - as.numeric(aux_cp_counts[6, ])
cp_counts <- cp_counts - as.numeric(aux_cp_counts[7, ])
# if the background subtraction is negative, assign the expression to 0
write.csv(cp_counts, quote = F, file = '/mnt/scratch2/STAMP/stamp11b_toRemoveBkg.csv')
# ran this on python and saved it as csv: /mnt/scratch2/STAMP/stamp11b_bkgRemoved.csv

no_bkg <- data.table::fread('/mnt/scratch2/STAMP/stamp11b_bkgRemoved.csv')
no_bkg <- as.data.frame(no_bkg)
no_bkg <- no_bkg[, -c(1:2)]
rownames(no_bkg) <- rownames(cp_counts)




centroids <- data.frame(x= cell_meta$y_global_px,
                        y = cell_meta$x_global_px,
                        cell = cell_meta$cell)
coords <- suppressWarnings(expr = CreateFOV(
  coords = centroids,
  type = 'centroids',
  key = 'FOV',
  assay = 'X-CP'
))

aux_assay <- CreateAssay5Object(counts = aux_cp_counts)

#stamp11b <- CreateSeuratObject(counts = Matrix::as.matrix(cp_counts, sparse = T), meta.data = cbind(cp_meta, cell_meta), assay = 'CP')
stamp11b <- CreateSeuratObject(counts = Matrix::as.matrix(no_bkg, sparse = T), meta.data = cbind(cp_meta, cell_meta), assay = 'CP')
stamp11b[['Controls']] <- aux_assay
suppressWarnings(expr = stamp11b[['fov']] <- coords) # add image to seurat object




## Finding out coordinate cuttoffs to separate the sub-stamps
ggplot(stamp11b[[]], aes(x_global_px,y_global_px)) + geom_point(size = 0.01) + coord_fixed()

# STAMP PBMC Mix is the larger square
stamp11b_pbmc <- subset(stamp11b,y_global_px > 40000 )
#

p <- ggplot(stamp11b_pbmc[[]], aes(x_global_px,y_global_px)) + geom_point(size = 0.001) + coord_fixed()
p + geom_hline(yintercept = 118000,  linetype = 'dashed', color = 'red', linewidth = 1.5) +
    geom_hline(yintercept = 52000,  linetype = 'dashed', color = 'red',linewidth = 1.5) +
    geom_vline(xintercept = 5000,  linetype = 'dashed', color = 'red',linewidth = 1.5) +
    geom_vline(xintercept = 71000,  linetype = 'dashed', color = 'red',linewidth = 1.5)

VlnPlot(stamp11b_pbmc, features = rownames(stamp11b_pbmc)[1:5], pt.size = 0, layer = 'counts')
VlnPlot(stamp11b_pbmc, features = c('nCount_CP', 'Area.um2', 'nFeature_CP'), pt.size = 0)
stamp11b_pbmc <- subset(stamp11b_pbmc, nCount_CP > 20 & nCount_CP < 2000 & Area.um2 < 250)





stamp11b_pbmc <- NormalizeData(stamp11b_pbmc, normalization.method = 'CLR', margin = 1)
stamp11b_pbmc <- ScaleData(stamp11b_pbmc, do.center = T)
stamp11b_pbmc <- FindVariableFeatures(stamp11b_pbmc)

stamp11b_pbmc <- RunPCA(stamp11b_pbmc)
stamp11b_pbmc <- FindNeighbors(stamp11b_pbmc, dims = 1:15)
stamp11b_pbmc <- FindClusters(stamp11b_pbmc, resolution = 0.5)
stamp11b_pbmc <- RunUMAP(stamp11b_pbmc, dims = 1:15)
DimPlot(stamp11b_pbmc, label = F, raster = F)
#FeaturePlot(stamp11b_pbmc, features = 'EpCAM', order = T, slot = 'scale.data') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#DotPlot(stamp11b_pbmc, features = rownames(stamp11b_pbmc), scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

#saveRDS(stamp11b_pbmc, file = 'stamp11b_pbmc_seurat.rds')

# Plot average expression per cluster
avg_stamp12b <- AverageExpression(stamp11b_pbmc, group.by = "seurat_clusters", layer = 'data')$CP
pheatmap(avg_stamp12b, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)
DotPlot(stamp11b_pbmc, features = rownames(avg_stamp12b), scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

# Lvl 1 from RNA (Ema): B    LowQ Myeloid      NK       T


stamp11b_pbmc$Lineage <- stamp11b_pbmc$seurat_clusters
stamp11b_pbmc[[]] <- stamp11b_pbmc[[]] %>%
  mutate(Lineage = case_when(Lineage %in% c('0','6', '9', '2', '5' ) ~ 'T',
                                   Lineage %in% c('10', '1', '8', '7') ~ 'Myeloid',
                                   Lineage %in% c('11') ~ 'LowQ',
                                   Lineage %in% c('3') ~ 'B',
                                   Lineage %in% c('4') ~ 'NK'))


png(width = 5, height = 5, filename = '/mnt/scratch2/STAMP/Figures/UMAP_Lineage_STAMP11b_XCP.png', units = 'in', res = 300)
DimPlot(stamp11b_pbmc, group.by = 'Lineage',label = F, raster = F, cols = lineage_colors) & NoAxes() & ggtitle('') & NoLegend()
dev.off()




#pdf("/mnt/scratch2/STAMP/Figures/STAMP11_BarPlot_ProportionClusterCellType_RNA_PRotein.pdf", 5, 10)
p6 <- ggplot(combined_df, aes(x = Dataset, y = Proportion, fill = CellType_pred)) +
  geom_bar(stat = "identity") +
  labs(x = "Dataset", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_bw() +theme(text=element_text(size=16)) +
  scale_y_continuous(labels = scales::percent_format())  + scale_fill_manual(values = brewer.pal(8, 'Paired')) +theme(axis.text.x = element_text(angle = 45, hjust = 1))




lineages <- unique(stamp11b_pbmc@meta.data$Lineage)
lineages <- lineages[1:4] # do not include lowq

# Create a list of Seurat objects, each subset by a category in Lineage
stamp11b_list <- lapply(lineages, function(lineage) {
  subset(stamp11b_pbmc, subset = Lineage == lineage)
})
names(stamp11b_list) <- lineages

stamp11b_list <- lapply(stamp11b_list, function(x) {
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 1)
  x <- ScaleData(x) %>% FindVariableFeatures()
  x <- RunPCA(x)
  x <- FindNeighbors(x, dims = 1:15)
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:15)
  return(x)
})

DimPlot(stamp11b_list$T)
DimPlot(stamp11b_list$B)
DimPlot(stamp11b_list$Myeloid)
DimPlot(stamp11b_list$NK)

## T

avg_stamp12b_T <- AverageExpression(stamp11b_list$T, group.by = "seurat_clusters", layer = 'data',
                                    features = c('4-1BB', 'CD127', 'CD27', 'CD3', "CD34", 'CD38', 'CD4',
                                                 'CD45RA', 'CD8','FOXP3', 'GITR', 'GZMA', 'GZMB', 'ICOS',
                                                 'LAG3', 'PD-1', 'TCF7', 'Tim-3'))$CP

pheatmap(avg_stamp12b_T, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)
DotPlot(stamp11b_list$T, features = rownames(stamp11b_list$T), scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

stamp11b_list$T$CellType <- stamp11b_list$T$seurat_clusters
stamp11b_list$T[[]] <- stamp11b_list$T[[]] %>%
  mutate(CellType = case_when(CellType %in% c('0', '5') ~ 'CD8 T',
                             CellType %in% c('0', '5', '4', '2', '7', '3', '1') ~ 'CD4 T',
                             CellType %in% c('6') ~ 'Naive T',
                             CellType %in% c('8', '10') ~ 'Activated T',
                             CellType %in% c('9') ~ 'Unresolved'
                             ))
t_cell_colors <- awtools::mpalette[1:5]
names(t_cell_colors) <- unique(stamp11b_list$T$CellType)

png(width = 5, height = 5, filename = '/mnt/scratch2/STAMP/Figures/UMAP_TCells_STAMP11b_XCP.png', units = 'in', res = 300)
DimPlot(stamp11b_list$T, group.by = 'CellType', raster = F, cols = t_cell_colors) & NoAxes() & ggtitle('') & NoLegend()
dev.off()


## Myeloid
DimPlot(stamp11b_list$Myeloid, label = T)

avg_stamp12b_Myeloid <-  AverageExpression(stamp11b_list$Myeloid, group.by = "seurat_clusters", layer = 'data')$CP
pheatmap(avg_stamp12b_Myeloid, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)

stamp11b_list$Myeloid$CellType <- stamp11b_list$Myeloid$seurat_clusters
stamp11b_list$Myeloid[[]] <- stamp11b_list$Myeloid[[]] %>%
  mutate(CellType = case_when(CellType %in% c('7') ~ 'DCs',
                              CellType %in% c('0','1','4', '6') ~ 'CD14 Monocytes',
                              CellType %in% c('3', '5', '8') ~ 'CD16 Monocytes',
                              CellType %in% c('2', '9') ~ 'Unresolved'))

DimPlot(stamp11b_list$Myeloid, group.by = 'CellType')

myeloid_cell_colors <- ggsci::pal_npg("nrc")(4)
names(myeloid_cell_colors) <- unique(stamp11b_list$Myeloid$CellType)

png(width = 5, height = 5, filename = '/mnt/scratch2/STAMP/Figures/UMAP_Myeloid_STAMP11b_XCP.png', units = 'in', res = 300)
DimPlot(stamp11b_list$Myeloid, group.by = 'CellType', raster = F, cols = myeloid_cell_colors) & NoAxes() & ggtitle('') & NoLegend()
dev.off()

# B
DimPlot(stamp11b_list$B, label = T)

avg_stamp12b_B <-  AverageExpression(stamp11b_list$B, group.by = "seurat_clusters", layer = 'data')$CP
pheatmap(avg_stamp12b_B, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)

stamp11b_list$B$CellType <- 'B'
png(width = 5, height = 5, filename = '/mnt/scratch2/STAMP/Figures/UMAP_B_STAMP11b_XCP.png', units = 'in', res = 300)
DimPlot(stamp11b_list$B, group.by = 'CellType', raster = F, cols = lineage_colors['B']) & NoAxes() & ggtitle('') & NoLegend()
dev.off()

# NK
DimPlot(stamp11b_list$NK, label = T)

avg_stamp12b_NK <-  AverageExpression(stamp11b_list$NK, group.by = "seurat_clusters", layer = 'data')$CP
pheatmap(avg_stamp12b_NK, scale = 'column', color = rev(brewer.pal(n = 11, name = "RdBu")), border_color = NA, cluster_rows = T, fontsize = 16, angle_col = 45)

stamp11b_list$NK$CellType <- 'NK'
png(width = 5, height = 5, filename = '/mnt/scratch2/STAMP/Figures/UMAP_NK_STAMP11b_XCP.png', units = 'in', res = 300)
DimPlot(stamp11b_list$NK, group.by = 'CellType', raster = F, cols = lineage_colors['NK']) & NoAxes() & ggtitle('') & NoLegend()
dev.off()



#
# ## Find CTCs based on Her2
# FeatureScatter(stamp11b_pbmc, feature1 = 'LAMP1', feature2 = 'EGFR', slot = 'data', raster = F, col = 'gray') + geom_hline(yintercept = 5, linetype ='dashed') +
#     geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs - MCF-7 in PBMCs')
#
# FeatureScatter(stamp11b_pbmc, feature1 = 'EGFR', feature2 = 'ICAM1', slot = 'counts', raster = F, col = 'gray') + geom_hline(yintercept = 5, linetype ='dashed') +
#   geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs - SK-BR-3 in PBMCs')
#
# # classify based on markers
# stamp11b_pbmc$CTC <- NA
# stamp11b_pbmc$CTC <- ifelse(
#     #stamp11b_pbmc@assays$CP$scale.data["EpCAM", ] > 2 &
#     stamp11b_pbmc@assays$CP$scale.data["EGFR", ] > 5 |
#     stamp11b_pbmc@assays$CP$scale.data["Her2", ] > 4,
#     #stamp11b_pbmc@assays$CP$scale.data["EpCAM", ] > 4,
#     'CTC',
#     'PBMC'
#   )
#
# # remove "CTCs" with small area (less than average PBMC cell)
# stamp11b_pbmc$CTC <- ifelse(
#   stamp11b_pbmc$CTC == "CTC" & stamp11b_pbmc$Area.um2 < median(stamp11b_pbmc$Area.um2),
#   "PBMC",  # Change label for "CTC" with area < mean
#   stamp11b_pbmc$CTC  # Keep the existing label if the condition is not met
# )
#
#
#
# VlnPlot(stamp11b_pbmc, features = c('Area.um2', 'nCount_CP'), group.by = 'CTC', pt.size = 0) & stat_compare_means()
# VlnPlot(stamp11b_pbmc, features = c('Her2', 'EGFR'), group.by = 'CTC', pt.size = 0, layer = 'scale.data') & stat_compare_means()
#
#
#
#
# FeatureScatter(stamp11b_pbmc, feature1 = 'Mean.PanCK', feature2 = 'EGFR', slot = 'counts', raster = F, group.by = 'CTC') + geom_hline(yintercept = 5, linetype ='dashed') +
#   geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs in PBMCs')
#
#
# stamp11b_pbmc[[]][stamp11b_pbmc$CTC %in% 'CTC', ]
#
# dittoBarPlot(stamp11b_pbmc, var = 'CTC', group.by = 'orig.ident')
#

##################################
sce_obj <- as.SingleCellExperiment(stamp11b_pbmc, assay = c("CP"))
writeH5AD(sce_obj, "/mnt/scratch2/STAMP/stamp11b_cp.h5ad", X_name = 'counts')






