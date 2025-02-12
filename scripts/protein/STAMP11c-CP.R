### STAMP 11c
library(tidyverse)
library(Seurat)
library(viridis)
library(plotly)

## Read in seurat manually

cp_counts <- data.table::fread('/mnt/scratch2/STAMP/stamp_11c/STAMPCP_exprMat_file.csv.gz')
cp_counts <- cp_counts[, -c(1:2)] # remove fov and cell id
aux_cp_counts <- cp_counts[, c(31:35,68:69)]
cp_counts <- cp_counts[, -c(31:35,68:69)] # remove membrane markers and control iGGs
cp_meta <- data.table::fread('/mnt/scratch2/STAMP/stamp_11c/STAMPCP_metadata_file.csv.gz')
cp_counts <- as.data.frame(t(cp_counts))
colnames(cp_counts) <- cp_meta$cell_id
cell_meta <- data.table::fread('/mnt/scratch2/STAMP/stamp_11c/STAMPCP-polygons.csv.gz') %>% as.data.frame
cell_meta$cell <- make.names(cell_meta$cell, unique = T)
cell_meta <- cell_meta[cell_meta$cell %in% colnames(cp_counts) , ]
aux_cp_counts <- as.data.frame(t(aux_cp_counts))
colnames(aux_cp_counts) <- cp_meta$cell_id

#Make data frame of only background
bkg_stamp11c <- data.frame('MsIgG1_STAMP11c' = t(aux_cp_counts[6,]),
                           'RbIgG_STAMP11c' = t(aux_cp_counts[7, ]))
bkg_stamp11c$STAMP <- 'STAMP11c'



## Remove background from counts (last 2 rows of aux_cp_counts)
#cp_counts <- cp_counts - as.numeric(aux_cp_counts[6, ])
#cp_counts <- cp_counts - as.numeric(aux_cp_counts[7, ])



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

stamp11c <- CreateSeuratObject(counts = Matrix::as.matrix(cp_counts, sparse = T), meta.data = cbind(cp_meta, cell_meta), assay = 'CP')
#stamp11b[['Controls']] <- aux_assay
#suppressWarnings(expr = stamp11b[['fov']] <- coords) # add image to seurat object




## Finding out coordinate cuttoffs to separate the sub-stamps
ggplot(stamp11c[[]], aes(x_global_px,y_global_px)) + geom_point(size = 0.01) + coord_fixed()

# STAMP PBMC Mix is the larger square
stamp11c_pbmc <- subset(stamp11c,y_global_px > 40000 )
#

p <- ggplot(stamp11c_pbmc[[]], aes(x_global_px,y_global_px)) + geom_point(size = 0.001) + coord_fixed()
p + geom_hline(yintercept = 118000,  linetype = 'dashed', color = 'red', linewidth = 1.5) +
  geom_hline(yintercept = 52000,  linetype = 'dashed', color = 'red',linewidth = 1.5) +
  geom_vline(xintercept = 5000,  linetype = 'dashed', color = 'red',linewidth = 1.5) +
  geom_vline(xintercept = 71000,  linetype = 'dashed', color = 'red',linewidth = 1.5)

VlnPlot(stamp11c_pbmc, features = rownames(stamp11c_pbmc)[1:5], pt.size = 0, layer = 'counts')
VlnPlot(stamp11c_pbmc, features = c('nCount_CP', 'Area.um2', 'nFeature_CP'), pt.size = 0)
#stamp11c_pbmc <- subset(stamp11c_pbmc, nCount_CP >= 0) ## 12.3% filtered (high background)





stamp11c_pbmc <- NormalizeData(stamp11c_pbmc)
stamp11c_pbmc <- ScaleData(stamp11c_pbmc)
stamp11c_pbmc <- FindVariableFeatures(stamp11c_pbmc)

stamp11c_pbmc <- RunPCA(stamp11c_pbmc)
stamp11c_pbmc <- FindNeighbors(stamp11c_pbmc, dims = 1:15)
stamp11c_pbmc <- FindClusters(stamp11c_pbmc, resolution = 0.5)
stamp11c_pbmc <- RunUMAP(stamp11c_pbmc, dims = 1:15)
DimPlot(stamp11c_pbmc)
FeaturePlot(stamp11c_pbmc, features = 'EpCAM', order = T, slot = 'scale.data') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
DotPlot(stamp11c_pbmc, features = rownames(stamp11c_pbmc), scale = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + RotatedAxis()

saveRDS(stamp11c_pbmc, file = 'stamp11c_pbmc_seurat.rds')


## Find CTCs based on Her2
FeatureScatter(stamp11b_pbmc, feature1 = 'LAMP1', feature2 = 'EGFR', slot = 'data', raster = F, col = 'gray') + geom_hline(yintercept = 5, linetype ='dashed') +
  geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs - MCF-7 in PBMCs')

FeatureScatter(stamp11b_pbmc, feature1 = 'EGFR', feature2 = 'ICAM1', slot = 'counts', raster = F, col = 'gray') + geom_hline(yintercept = 5, linetype ='dashed') +
  geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs - SK-BR-3 in PBMCs')

# classify based on markers
stamp11b_pbmc$CTC <- NA
stamp11b_pbmc$CTC <- ifelse(
  #stamp11b_pbmc@assays$CP$scale.data["EpCAM", ] > 2 &
  stamp11b_pbmc@assays$CP$scale.data["EGFR", ] > 5 |
    stamp11b_pbmc@assays$CP$scale.data["Her2", ] > 4,
  #stamp11b_pbmc@assays$CP$scale.data["EpCAM", ] > 4,
  'CTC',
  'PBMC'
)

# remove "CTCs" with small area (less than average PBMC cell)
stamp11b_pbmc$CTC <- ifelse(
  stamp11b_pbmc$CTC == "CTC" & stamp11b_pbmc$Area.um2 < median(stamp11b_pbmc$Area.um2),
  "PBMC",  # Change label for "CTC" with area < mean
  stamp11b_pbmc$CTC  # Keep the existing label if the condition is not met
)



VlnPlot(stamp11b_pbmc, features = c('Area.um2', 'nCount_CP'), group.by = 'CTC', pt.size = 0) & stat_compare_means()
VlnPlot(stamp11b_pbmc, features = c('Her2', 'EGFR'), group.by = 'CTC', pt.size = 0, layer = 'scale.data') & stat_compare_means()




FeatureScatter(stamp11b_pbmc, feature1 = 'Mean.PanCK', feature2 = 'EGFR', slot = 'counts', raster = F, group.by = 'CTC') + geom_hline(yintercept = 5, linetype ='dashed') +
  geom_vline(xintercept = 70, linetype ='dashed') + ggtitle('CTCs in PBMCs')


stamp11b_pbmc[[]][stamp11b_pbmc$CTC %in% 'CTC', ]

dittoBarPlot(stamp11b_pbmc, var = 'CTC', group.by = 'orig.ident')


##################################

