# Save plots at /mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Figure7

setwd("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Figure7")
figures_dir <- "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Figure7"
## Load packatges
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)

### Load data
Sobj_merg_CellNuc <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_merg_CellNuc_HiLow_cells.rds")

Sobj_merg_CellNuc@meta.data <- Sobj_merg_CellNuc@meta.data %>%
  mutate(Tissue_type = case_when(
    grepl("BRAIN|brain|HEART|heart|KIDNEY|kidney|LIVER|liver|LUNG|lung|MIX|mix", Well) ~
      ifelse(grepl("BRAIN|brain", Well), "brain",
             ifelse(grepl("HEART|heart", Well), "heart",
                    ifelse(grepl("KIDNEY|kidney", Well), "kidney",
                           ifelse(grepl("LIVER|liver", Well), "liver",
                                  ifelse(grepl("LUNG|lung", Well), "lung", ifelse(grepl("MIX|mix", Well), "mix", 
                                                                                  NA)))))),
    TRUE ~ NA_character_ # Default case for other tissues
  ))


### Define functions
Filter_Query_SCdata <- function(Sobj = NULL,
                                Organ = NULL,
                                Cell_Nuc = NULL,
                                High_Low_QC = NULL) {
  print("Filtering by Organ")
  # Check if Organ is NULL
  if (!is.null(Organ)) {
    Sobj <- subset(Sobj, Tissue_type %in% Organ)
  }
  
  print("Filtering by Cell_Nuc")
  # Check if Cell_Nuc is NULL
  if (!is.null(Cell_Nuc)) {
    Sobj <- subset(Sobj, Cell_or_Nuc %in% Cell_Nuc)
  }
  
  print("Filtering by High_Low_QC")
  # Check if High_Low_QC is NULL
  if (!is.null(High_Low_QC)) {
    Sobj <- subset(Sobj, High_low_QCcells %in% High_Low_QC)
  }
  
  # Join assay layers
  Sobj@assays$RNA <- Sobj@assays$Xenium
  Sobj[["RNA"]] <- JoinLayers(Sobj[["RNA"]])
  
  print(paste0("Query data has dim of: ", dim(Sobj)))
  return(Sobj)
}


### Figure 7  ---------------------------------------------------- 

### QC Barplots & Scatter Plot ----------------

# Proportion Barplot - high low QC cells/nucs ------------
metadata <- Sobj_merg_CellNuc@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality cells"
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality cells"

metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality nucs"

table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_cells_and_nucs.csv")

png(width = 5, height = 5, filename = paste0(figures_dir,'/CellProp_QC_Barplot_NoLegend.png'), units = 'in', res = 300)


plot <- ggplot(barplot_inp, 
               aes(x = Cell_or_Nuc, y = proportion, fill = High_low_QCcells_2)) + 
  geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
  scale_fill_manual(values = c("High quality cells" =  "darkblue", 
                               "Low quality cells" = "lightblue3",
                               "High quality nucs" = "darkred",
                               "Low quality nucs" = "darksalmon")) +
  # scale_color_manual(values = c("brain" = "blue", 
  #                               "heart" = "orange", 
  #                               "kidney" = "purple",
  #                               "liver" = "black",
  #                               "lung" = "brown",
  #                               "mix" = "gray")) + # Adjust colors as needed
  labs(title = "",
       x = "",
       y = "Proportion of Cells (%)") + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  facet_grid(~ Tissue_type) + NoLegend()

# Display the plot
print(plot)

dev.off()


# Scatter plot -------------
Sobj_cells <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = NULL,
                                  Cell_Nuc = "Cells",
                                  High_Low_QC = "High quality cells")
Sobj_cells <- NormalizeData(Sobj_cells)
cells_exp <- FetchData(Sobj_cells, vars = rownames(Sobj_cells))
cells_exp[1:4, 1:4]

genes_cells_mean_exp <- colMeans(cells_exp, na.rm = TRUE) %>% data.frame()
genes_cells_mean_exp %>% head()
dim(genes_cells_mean_exp)
names(genes_cells_mean_exp) <- "cells_mean_gene_exp"

Sobj_nucs <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                 Organ = NULL,
                                 Cell_Nuc = "Nucs",
                                 High_Low_QC = "High quality cells")
Sobj_nucs <- NormalizeData(Sobj_nucs)
nucs_exp <- FetchData(Sobj_nucs, vars = rownames(Sobj_nucs))
nucs_exp[1:4, 1:4]

genes_nucs_mean_exp <- colMeans(nucs_exp, na.rm = TRUE) %>% data.frame()
genes_nucs_mean_exp %>% head()
dim(genes_nucs_mean_exp)
names(genes_nucs_mean_exp) <- "nucs_mean_gene_exp"

identical(rownames(genes_nucs_mean_exp),
          rownames(genes_cells_mean_exp))

mean_gene_exp <- cbind(genes_cells_mean_exp,genes_nucs_mean_exp)
# saveRDS(mean_gene_exp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/AvgExpCellsVsNucs_df.rds")

mean_gene_exp <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/AvgExpCellsVsNucs_df.rds")

png(width = 5, height = 5, filename = paste0(figures_dir,'/Correlation_MeanGeneExp_cells_vs_nucs_scatterplot.png'), units = 'in', res = 300)
library(ggplot2); #theme_set(theme_classic())
library(ggpubr)

ggplot(mean_gene_exp, aes(x = cells_mean_gene_exp, y = nucs_mean_gene_exp)) + 
  geom_point(aes(size = I(3), alpha = 0.5), shape = 16, color = "gray") + 
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") + # Add trend line
  labs(title = "",
       x = "Average Expression Cells",
       y = "Average Expression Nuclei") +
  stat_cor(method = "pearson", label.x = 0, label.y = 4) +
  theme_minimal() +
  #scale_y_break(c(3000, 8755)) +
  #ylim(0, 8755) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), 
    strip.text.x = element_text(size = 16)
  ) + NoLegend()
  

dev.off()


# 3 Variable QC boxplot ---------------

png(width = 5, height = 5, filename = paste0(figures_dir,'/nCount_Boxplot_NoLegend.png'), units = 'in', res = 300)

# library(ggplot2)
# install.packages("ggbreak")
# library(ggbreak)

library(ggplot2)
ggplot(metadata, aes(x = Cell_or_Nuc, y = log(nCount_Xenium + 1), fill = High_low_QCcells_2)) +
  geom_boxplot(outlier.color = "white"
    #outlier.shape = NA
    ) +
  scale_fill_manual(values = c("High quality cells" =  "darkblue", 
                               "Low quality cells" = "lightblue3",
                               "High quality nucs" = "darkred",
                               "Low quality nucs" = "darksalmon")) +
  labs(title = "",
       x = "",
       y = "nCount") +
  theme_minimal() +
  #scale_y_break(c(3000, 8755)) +
  #ylim(0, 8755) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), 
    strip.text.x = element_text(size = 16)
  ) + facet_grid(~ Tissue_type, scales = "free", ) + NoLegend()

dev.off()


png(width = 5, height = 5, filename = paste0(figures_dir,'/nFeature_Boxplot_NoLegend.png'), units = 'in', res = 300)

library(ggplot2)
ggplot(metadata, aes(x = Cell_or_Nuc, y = log(nFeature_Xenium + 1), fill = High_low_QCcells_2)) +
  geom_boxplot(outlier.color = "white"
               #outlier.shape = NA
  ) +
  scale_fill_manual(values = c("High quality cells" =  "darkblue", 
                               "Low quality cells" = "lightblue3",
                               "High quality nucs" = "darkred",
                               "Low quality nucs" = "darksalmon")) +
  labs(title = "",
       x = "",
       y = "nFeature") +
  theme_minimal() +
  #scale_y_break(c(3000, 8755)) +
  #ylim(0, 8755) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), 
    strip.text.x = element_text(size = 16)
  ) + facet_grid(~ Tissue_type, scales = "free", ) + NoLegend()

dev.off()



png(width = 5, height = 5, filename = paste0(figures_dir,'/CellArea_Boxplot_NoLegend.png'), units = 'in', res = 300)

library(ggplot2)
ggplot(metadata, aes(x = Cell_or_Nuc, y = log(cell_area + 1), fill = High_low_QCcells_2)) +
  geom_boxplot(outlier.color = "white"
               #outlier.shape = NA
  ) +
  scale_fill_manual(values = c("High quality cells" =  "darkblue", 
                               "Low quality cells" = "lightblue3",
                               "High quality nucs" = "darkred",
                               "Low quality nucs" = "darksalmon")) +
  labs(title = "",
       x = "",
       y = "Cell Area") +
  theme_minimal() +
  #scale_y_break(c(3000, 8755)) +
  #ylim(0, 8755) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), 
    strip.text.x = element_text(size = 16)
  ) + facet_grid(~ Tissue_type, scales = "free", ) + NoLegend()

dev.off()


### Cells UMAPs 

## Individual Organs - Cell Type UMAP
Sobj_organs_annotted_NOmix <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")
# Define vectors for organs and modalities
Organs <- c("Heart", "Lung", "Liver", "Kidney", "Brain")  # Add more organs as needed
Modalities <- c("Cells", "Nuclei")    # Add more modalities as needed

Modality <- Modalities[1] # keep it 1 for "Cells"
Organ <- Organs[5] # PLOT it from 1 to 5

# Open a PNG device
png(width = 5, height = 5, filename = paste0(figures_dir,'/', Organ, "_", Modality ,'_Celltype_UMAP_NoLegend.png'), units = 'in', res = 300)

# Create the UMAP plot with error handling
celltypes <- Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query$predicted.id %>% table() %>% names()
DimPlot(Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query, 
        reduction = "ref.umap", 
        group.by = 'predicted.id', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(celltypes))) & NoAxes() & ggtitle('') & NoLegend()

# Close the PNG device
dev.off()


png(width = 5, height = 5, filename = paste0(figures_dir,'/', Organ, "_", Modality ,'_Celltype_UMAP.png'), units = 'in', res = 300)

DimPlot(Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query, 
        reduction = "ref.umap", 
        group.by = 'predicted.id', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(celltypes))) & NoAxes() & ggtitle('') 


dev.off()


## Cells Mix - celltypes  
data.query <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Cells_Mix_Predicted_Orngan_Celltype.rds")

png(width = 5, height = 5, filename = paste0(figures_dir,'/', 'Mix', "_", 'Cells' ,'_Celltype_UMAP.png'), units = 'in', res = 300)

# celltypes <- data.query$predicted.celltype %>% table() %>% names()
# Collor pallet for > 10 categories 
tableau_colors <- tableau_color_pal("Tableau 20")(20)
library(RColorBrewer)
additional_colors <- brewer.pal(n = 10, name = "Set3")
combined_colors <- c(tableau_colors, additional_colors)

DimPlot(data.query, 
        reduction = "ref.umap", 
        group.by = 'predicted.celltype', 
        raster = F, 
        cols = combined_colors) & NoAxes() & ggtitle('') #& NoLegend()

dev.off()


## Cells Mix - Organs  
png(width = 5, height = 5, filename = paste0(figures_dir,'/', 'Mix', "_", 'Cells' ,'_Organ_UMAP_NoLegend.png'), units = 'in', res = 300)

organs_name <- data.query$predicted.Well %>% table() %>% names()
DimPlot(data.query, 
        reduction = "ref.umap", 
        group.by = 'predicted.Well', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(organs_name))) & NoAxes() & ggtitle('') & NoLegend()

dev.off()


## Cell Mix - Louvain Cluster ~ Organ

# # Collor pallet for > 10 categories 
# combined_colors <- c(
#   "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#   "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
#   "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
#   "#ccebc5", "#ffed6f", "#ff7f00", "#c2c2f0", "#ffb3e6",
#   "#ffcc99", "#66b3ff", "#99ff99", "#ffcc00",
#   "#cce5ff", "#ffe6cc", "#cc0033",
#   "#999999", "#663399",
#   "#ff6699",
#   "#999966",
#   "#ff9966"
# )
# 
# # Remove duplicates
# unique_colors <- unique(combined_colors)
# 
# # Add two new colors
# new_colors <- c("#ffcc33", "#33ccff")
# combined_colors <- c(unique_colors, new_colors)

png(width = 3, height = 7, filename = paste0(figures_dir,'/', 'Mix', "_", 'Cells' ,'_Mix_vs_Organ_PropBarplot_NoLegend.png'), units = 'in', res = 300)

organs_name <- data.query$predicted.Well %>% table() %>% names()
organ_colors <- ggthemes::tableau_color_pal()(length(organs_name))
library(dittoSeq)
dittoBarPlot(data.query, 
             "predicted.Well", 
             group.by = "Tissue_type", 
             colors = seq_along(organ_colors),
             color.panel = organ_colors) + 
  labs(title = "",
       x = "Organ in Mix cells",
       y = "Percent of Cells") +
  scale_x_discrete(labels = c("mix" = "")) +
  theme(axis.text.x = element_text(angle = 0, hjust = -1, size=16),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) & ggtitle('')  & NoLegend()

dev.off()



## Dotplot -  Organ markers 
#  Define cluster markers
library(presto)
markers <- wilcoxauc(data.query, 'predicted.Well')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 6) %>% 
  ungroup() -> top_markers
table(top_markers$group)
df <- data.frame(top_markers)

data.query$predicted.Well <- factor(data.query$predicted.Well, levels =c("BRAIN_cells", "HEART_cells","KIDNEY_cells", "LIVER_cells","LUNG_cells")) #reordenando os leveis

library(plyr)
data.query$predicted.Well <- mapvalues(data.query$predicted.Well, from = c("BRAIN_cells", "HEART_cells","KIDNEY_cells", "LIVER_cells","LUNG_cells"), to = c("Brain", "Heart","Kidney", "Liver","Lung")) #alterando o nome dos leveis


png(width = 12, height = 5, filename = paste0(figures_dir,'/', 'Mix', "_", 'Cells' ,'_Organs_Markers_Dotplot.png'), units = 'in', res = 300)

# Plot Dotplot 
library(ggplot2)
library(viridis)
Idents(data.query) <- "predicted.Well"
DotPlot(data.query, 
        features =  unique(c(df$feature))
        #features =  unique(c(top_markers_auc))
) +
  labs(title = "", x = "Organ Gene Markers", y = "") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) + ggtitle("") + NoLegend()

dev.off()

















### Supp figures ----------------------------------------------------
figures_dir <- "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Figure7/SuppFig"
## embedings (both cells and nucs) ------------
## ncount vs nfeature full dataset ~ facet -----------

Sobj_hiqc <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                 Organ = NULL,
                                 Cell_Nuc = NULL,
                                 High_Low_QC = "High quality cells")

png(width = 10, height = 5, filename = paste0(figures_dir,'/', 'OrgansMix', "_", 'CellsNucs' ,'_highQC_nCounter_vs_nFeature_scatter.png'), units = 'in', res = 300)


metadata <- Sobj_hiqc@meta.data
library(ggplot2); theme_set(theme_classic())
plt <- ggplot(metadata, aes(x = nCount_Xenium, 
                     y = nFeature_Xenium, 
                     #shape = Cell_or_Nuc,
                     color = Cell_or_Nuc)) +
  scale_color_manual(values = c("Cells" = "darkred", 
                                "Nucs" = "lightblue")) +  # Closing parenthesis added here
  geom_point(size = 2, alpha = 0.5) + 
  labs(title = "",
       x = "nCount",
       y = "nFeature") +
  theme_minimal() +
  #scale_y_break(c(3000, 8755)) +
  #ylim(0, 8755) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), 
    strip.text.x = element_text(size = 18)
  ) + facet_grid(~ Tissue_type) & ggtitle('')  #& NoLegend()

print(plt)

dev.off()


## Nucs UMAPs ------------
## Individual Organs - Cell Type UMAP
Sobj_organs_annotted_NOmix <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")
# Define vectors for organs and modalities
Organs <- c("Heart", "Lung", "Liver", "Kidney", "Brain")  # Add more organs as needed
Modalities <- c("Cells", "Nuclei")    # Add more modalities as needed

Modality <- Modalities[2] # keep it 1 for "Cells"
Organ <- Organs[5] # PLOT it from 1 to 5

# Open a PNG device
png(width = 5, height = 5, filename = paste0(figures_dir,'/', Organ, "_", Modality ,'_Celltype_UMAP_NoLegend.png'), units = 'in', res = 300)

# Create the UMAP plot with error handling
celltypes <- Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query$predicted.id %>% table() %>% names()
DimPlot(Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query, 
        reduction = "ref.umap", 
        group.by = 'predicted.id', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(celltypes))) & theme(text=element_text(size=16)) & NoAxes() & ggtitle('') & NoLegend()

# Close the PNG device
dev.off()


png(width = 5, height = 5, filename = paste0(figures_dir,'/', Organ, "_", Modality ,'_Celltype_UMAP.png'), units = 'in', res = 300)

DimPlot(Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query, 
        reduction = "ref.umap", 
        group.by = 'predicted.id', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(celltypes))) & theme(text=element_text(size=16)) & NoAxes() & ggtitle('') 


dev.off()




## Nucs Mix - celltypes 


## Merge celltypes 
Sobj_Nucs_Mix_annotted <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_Nucs_Mix_annotted_NOorgans.rds")
Sobj_Nucs_Mix_annotted_2 <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_Nucs_Mix_annotted_NOorgans_2.rds") # it has the cell type projected on it from another run I did on terminal

length(intersect(colnames(Sobj_Nucs_Mix_annotted$Sobj_query),
                 colnames(Sobj_Nucs_Mix_annotted_2$Sobj_query)))

Sobj_Nucs_Mix_annotted <- Sobj_Nucs_Mix_annotted$Sobj_query
Sobj_Nucs_Mix_annotted_2 <- Sobj_Nucs_Mix_annotted_2$Sobj_query

names(Sobj_Nucs_Mix_annotted_2@meta.data)[names(Sobj_Nucs_Mix_annotted_2@meta.data) == "predicted.id"] <- "predicted.celltype"
names(Sobj_Nucs_Mix_annotted_2@meta.data)[names(Sobj_Nucs_Mix_annotted_2@meta.data) == "predicted.id.score"] <- "predicted.celltype.score"


names(Sobj_Nucs_Mix_annotted@meta.data)[names(Sobj_Nucs_Mix_annotted@meta.data) == "predicted.id"] <- "predicted.Well"
names(Sobj_Nucs_Mix_annotted@meta.data)[names(Sobj_Nucs_Mix_annotted@meta.data) == "predicted.id.score"] <- "predicted.Well.score"



df_meta <- Sobj_Nucs_Mix_annotted_2@meta.data
df_meta <- df_meta[, c("predicted.celltype", "predicted.celltype.score")]
data.query = AddMetaData(
  object = Sobj_Nucs_Mix_annotted,
  metadata = df_meta)
data.query_nuc <- data.query
# saveRDS(data.query_nuc, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Nucs_Mix_Predicted_Orngan_Celltype.rds")
data.query <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Nucs_Mix_Predicted_Orngan_Celltype.rds")


png(width = 6, height = 6, filename = paste0(figures_dir,'/', 'Mix', "_", 'Nucs' ,'_Celltype_UMAP.png'), units = 'in', res = 300)

# celltypes <- data.query$predicted.celltype %>% table() %>% names()
# Collor pallet for > 10 categories 
library(ggthemes)
tableau_colors <- tableau_color_pal("Tableau 20")(20)
library(RColorBrewer)
additional_colors <- brewer.pal(n = 10, name = "Set3")
combined_colors <- c(tableau_colors, additional_colors)

DimPlot(data.query, 
        reduction = "ref.umap", 
        group.by = 'predicted.celltype', 
        raster = F, 
        cols = combined_colors) & NoAxes() & ggtitle('') #& NoLegend()

dev.off()


## Cells Mix - Organs  
png(width = 5, height = 5, filename = paste0(figures_dir,'/', 'Mix', "_", 'Nucs' ,'_Organ_UMAP.png'), units = 'in', res = 300)

organs_name <- data.query$predicted.Well %>% table() %>% names()
DimPlot(data.query, 
        reduction = "ref.umap", 
        group.by = 'predicted.Well', 
        raster = F, 
        cols = ggthemes::tableau_color_pal()(length(organs_name))) & NoAxes() & ggtitle('') #& NoLegend()

dev.off()


## Cell Mix - Louvain Cluster ~ Organ

# # Collor pallet for > 10 categories 
# combined_colors <- c(
#   "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#   "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
#   "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
#   "#ccebc5", "#ffed6f", "#ff7f00", "#c2c2f0", "#ffb3e6",
#   "#ffcc99", "#66b3ff", "#99ff99", "#ffcc00",
#   "#cce5ff", "#ffe6cc", "#cc0033",
#   "#999999", "#663399",
#   "#ff6699",
#   "#999966",
#   "#ff9966"
# )
# 
# # Remove duplicates
# unique_colors <- unique(combined_colors)
# 
# # Add two new colors
# new_colors <- c("#ffcc33", "#33ccff")
# combined_colors <- c(unique_colors, new_colors)

png(width = 3, height = 7, filename = paste0(figures_dir,'/', 'Mix', "_", 'Nucs' ,'_Mix_vs_Organ_PropBarplot_NoLegend.png'), units = 'in', res = 300)

organ_colors <- ggthemes::tableau_color_pal()(length(organs_name))
library(dittoSeq)
dittoBarPlot(data.query, 
             "predicted.Well", 
             group.by = "Tissue_type", 
             colors = seq_along(organ_colors),
             color.panel = organ_colors) + 
  labs(title = "",
       x = "Organ in Mix",
       y = "Percent of Cells") +
  scale_x_discrete(labels = c("mix" = "")) +
  theme(axis.text.x = element_text(angle = 0, hjust = -1, size=16),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) & ggtitle('')  & NoLegend()

dev.off()





## Dotplot -  Organ markers 
#  Define cluster markers
library(presto)
markers <- wilcoxauc(data.query, 'predicted.Well')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 6) %>% 
  ungroup() -> top_markers
table(top_markers$group)
df <- data.frame(top_markers)

data.query$predicted.Well <- factor(data.query$predicted.Well, levels =c("brain_nuc", "heart_nuc","kidney_nuc", "liver_nuc","lung_nuc")) #reordenando os leveis
library(plyr)
data.query$predicted.Well <- mapvalues(data.query$predicted.Well, from = c("brain_nuc", "heart_nuc","kidney_nuc", "liver_nuc","lung_nuc"), to = c("Brain", "Heart","Kidney", "Liver","Lung"))

png(width = 12, height = 5, filename = paste0(figures_dir,'/', 'Mix', "_", 'Nucs' ,'_Organs_Markers_Dotplot_NoLegend.png'), units = 'in', res = 300)

# Plot Dotplot 
library(ggplot2)
library(viridis)
Idents(data.query) <- "predicted.Well"
DotPlot(data.query, 
        features =  unique(c(df$feature))
        #features =  unique(c(top_markers_auc))
) +
  labs(title = "", x = "Organ Gene Markers", y = "") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) + ggtitle("") + NoLegend()

dev.off()















### Tables ---------

## All cells and nucs (no cell type in here)
metadata <- Sobj_merg_CellNuc@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality cells"
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality cells"

metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality nucs"

table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_ALL_cells_and_nucs.csv")


## Organs - Cells - cell types 
Sobj_organs_annotted_NOmix <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")

obj1 <- Sobj_organs_annotted_NOmix[["Cells"]][["Heart"]]$Sobj_query
obj2 <- Sobj_organs_annotted_NOmix[["Cells"]][["Lung"]]$Sobj_query
obj3 <- Sobj_organs_annotted_NOmix[["Cells"]][["Liver"]]$Sobj_query
obj4 <- Sobj_organs_annotted_NOmix[["Cells"]][["Kidney"]]$Sobj_query
obj5 <- Sobj_organs_annotted_NOmix[["Cells"]][["Brain"]]$Sobj_query

obj_cells <- Reduce(function(x, y) merge(x, y), list(obj1,
                                                     obj2,
                                                     obj3,
                                                     obj4,
                                                     obj5))

names(obj_cells@meta.data)[names(obj_cells@meta.data) == "predicted.id"] <- "predicted.celltype"

metadata <- obj_cells@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"


table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type, predicted.celltype) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_organs_cells_celltype.csv")


## Mix - Cells - cell types
data.query <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Cells_Mix_Predicted_Orngan_Celltype.rds")
metadata <- data.query@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality cells"
# metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality cells"

# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality nucs"

table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type, predicted.celltype) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_mix_cells_celltype.csv")



## Organs - Nucs - cell types 
Sobj_organs_annotted_NOmix <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")

obj1 <- Sobj_organs_annotted_NOmix[["Nuclei"]][["Heart"]]$Sobj_query
obj2 <- Sobj_organs_annotted_NOmix[["Nuclei"]][["Lung"]]$Sobj_query
obj3 <- Sobj_organs_annotted_NOmix[["Nuclei"]][["Liver"]]$Sobj_query
obj4 <- Sobj_organs_annotted_NOmix[["Nuclei"]][["Kidney"]]$Sobj_query
obj5 <- Sobj_organs_annotted_NOmix[["Nuclei"]][["Brain"]]$Sobj_query

obj_nucs <- Reduce(function(x, y) merge(x, y), list(obj1,
                                                    obj2,
                                                    obj3,
                                                    obj4,
                                                    obj5))
names(obj_nucs@meta.data)[names(obj_nucs@meta.data) == "predicted.id"] <- "predicted.celltype"

metadata <- obj_nucs@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
# metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality cells"

# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality nucs"

table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type, predicted.celltype) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_organs_nucs_celltype.csv")



## Mix - Nucs - cell types
data.query <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Nucs_Mix_Predicted_Orngan_Celltype.rds")

metadata <- data.query@meta.data
metadata$High_low_QCcells_2 <- NA
metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
# metadata[metadata$Cell_or_Nuc %in% "Cells" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality cells"

# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "High quality cells", ]$High_low_QCcells_2 <- "High quality nucs"
# metadata[metadata$Cell_or_Nuc %in% "Nucs" & metadata$High_low_QCcells %in% "Low quality cells", ]$High_low_QCcells_2 <- "Low quality nucs"

table(metadata$High_low_QCcells_2)

library(dplyr)
quality_counts <- metadata %>%
  group_by(High_low_QCcells_2, Well, Cell_or_Nuc, Tissue_type, predicted.celltype) %>%  # Group by both High_low_QCcells and Well
  dplyr::summarise(count = n(), .groups = 'drop') %>%  # Count occurrences and drop grouping
  group_by(Well) %>% data.frame()

out_list <- list()
for(Well in levels(as.factor(quality_counts$Well))) {
  barplot_inp <- quality_counts[quality_counts$Well %in% Well, ] %>% 
    mutate(
      total_count = sum(count),  # Calculate total count per Well (only for the same Well)
      proportion = (count / total_count) * 100)
  
  out_list[[Well]] <- barplot_inp
}

barplot_inp <- do.call(rbind, out_list)

# write.csv(barplot_inp, "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_mix_nucs_celltype.csv")





### Number of cells on Reference scRNA-seq 
setwd("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS")

brain_1 <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Brain_Microglia-counts.csv") %>% data.frame()
brain_2 <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Brain_Neurons-counts.csv") %>% data.frame()
dim(brain_2)[2] + dim(brain_1)[2]


heart <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Heart-counts.csv") %>% data.frame()

liver <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Liver-counts.csv") %>% data.frame()


lung <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Lung-counts.csv") %>% data.frame()

kidney <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/Kidney-counts.csv") %>% data.frame()




data.frame(organ = c("brain", "heart", "liver", "lung", "kidney"),
           nCells = c(c(dim(brain_2)[2] + dim(brain_1)[2]),
                      dim(heart)[2],
                      dim(liver)[2],
                      dim(lung)[2],
                      dim(kidney)[2]))



# Load necessary libraries
library(readr)
library(dplyr)

# Define the file names
file_names <- c("Bladder-counts.csv", "Brain_Microglia-counts.csv", 
                "Brain_Neurons-counts.csv", "Colon-counts.csv", 
                "Fat-counts.csv", "Heart-counts.csv", 
                "Kidney-counts.csv", "Liver-counts.csv", 
                "Lung-counts.csv", "Mammary-counts.csv", 
                "Marrow-counts.csv", "Muscle-counts.csv", 
                "Pancreas-counts.csv", "Skin-counts.csv", 
                "Spleen-counts.csv", "Thymus-counts.csv", 
                "Tongue-counts.csv", "Trachea-counts.csv")

# Initialize an empty data frame to store the results
column_counts <- data.frame(File = character(), Columns = integer(), stringsAsFactors = FALSE)

# Loop through each file, read it, and get the number of columns
for (file in file_names) {
  # Construct the full path to the file (update this path as necessary)
  file_path <- paste0("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS/", file)
  
  # Read the CSV file into a data frame
  df <- readr::read_csv(file_path) %>% data.frame()
  
  # Get the number of columns and append to the results data frame
  column_counts <- rbind(column_counts, data.frame(File = file, Columns = ncol(df)))
}

# Print the results
print(column_counts)
column_counts$Columns %>% sum()



library(readr)
NumberOfCells_ALL_cells_and_nucs <- readr::read_csv("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/Tables/NumberOfCells_ALL_cells_and_nucs.csv")
df <- data.frame(NumberOfCells_ALL_cells_and_nucs)


organs <- c("brain", "heart", "kidney", "liver", "lung", "mix")
# Initialize a data frame to store counts
counts_df <- data.frame(Organ = organs, Count = integer(length(organs)))
# Loop through each organ to calculate counts
for (i in seq_along(organs)) {
  counts_df$Count[i] <- df[df$Tissue_type %in% organs[i] & df$High_low_QCcells_2 %in% "High quality cells", ]$count
}
# Print the resulting counts table
counts_df_cells <- counts_df
counts_df_cells$Prep <- "Cell"


# Initialize a data frame to store counts
counts_df <- data.frame(Organ = organs, Count = integer(length(organs)))
# Loop through each organ to calculate counts
for (i in seq_along(organs)) {
  counts_df$Count[i] <- df[df$Tissue_type %in% organs[i] & df$High_low_QCcells_2 %in% "High quality nucs", ]$count
}
# Print the resulting counts table
counts_df_nucs <- counts_df
counts_df_nucs$Prep <- "Nucs"

print(counts_df_cells)
print(counts_df_nucs)


