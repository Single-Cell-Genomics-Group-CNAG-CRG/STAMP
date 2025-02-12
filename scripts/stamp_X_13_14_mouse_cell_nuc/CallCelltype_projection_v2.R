# Cell type reference from 
# https://github.com/czbiohub-sf/tabula-muris (github)
# https://pubmed.ncbi.nlm.nih.gov/30283141/ (paper)
# https://explore.data.humancellatlas.org/projects/e0009214-c0a0-4a7b-96e2-d6a83e966ce0 (single cell portal)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109774 (GSE)

# Description: Using FACs single cell data as cell type reference to annotate STAMP_Xenium_Mouse data

# Cell type ref data at
# setwd("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS")


### Load pacakges 
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

ref_data_dir <- "/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/FACS"

metadata <- readr::read_csv("/mnt/scratch2/Maycon/STAMPX_13_14_mouse_5k/Mouse_atlas_sccelltype/FACS/annotations_FACS.csv") %>% data.frame()
rownames(metadata) <- metadata$cell

### Define functions 

Read_Ref_SCdata <- function(data_dir = NULL, 
                            cell_type_pattern = NULL, 
                            metadata = NULL) {
  # List all subdirectories in the specified data directory
  all_files <- list.files(data_dir, full.names = TRUE)
  
  # Initialize an empty list to store Seurat objects for each matching directory
  seurat_objects <- list()
  
  for(i in 1:length(cell_type_pattern)) {
  # Loop through each directory and read data if it matches the cell type pattern
  if (any(grepl(cell_type_pattern[i], basename(all_files)))) {
    
    mtx_path <- all_files[grepl(cell_type_pattern[i], all_files)]
    mtx <- data.table::fread(mtx_path) %>% data.frame()
    rownames(mtx) <- mtx$V1
    mtx$V1 <- NULL
    
    # Create a Seurat object for this dataset
    seurat_object <- CreateSeuratObject(counts = mtx, 
                                        project = cell_type_pattern[i])
    
  }
    
    # Subset metadata by cell type 
    metadata <- metadata[metadata$tissue %in% cell_type_pattern[i], ]
    
    # Check if row names of metadata match with cell barcodes
    if (any(rownames(metadata) %in% colnames(seurat_object))) {
      
      # Subset metadata to match only the cells present in this Seurat object
      seurat_object <- AddMetaData(seurat_object, 
                                   metadata[rownames(metadata) %in% colnames(seurat_object), ])
    }
    
    remove_cells <- seurat_object@meta.data[seurat_object@meta.data$cell_ontology_class %in% NA | seurat_object@meta.data$cell_ontology_class %in% "unknown", ] %>% rownames()
    seurat_object <- subset(x = seurat_object, cells = remove_cells, invert = TRUE)

    
    # Append the Seurat object to the list
    seurat_objects[[cell_type_pattern[i]]] <- seurat_object
    
  }
    # Combine all Seurat objects into one
    combined_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)
    combined_seurat_object[["RNA"]] <- JoinLayers(combined_seurat_object[["RNA"]])
    
    print(paste0("Ref data has dim of: ", dim(combined_seurat_object)))
    return(combined_seurat_object)
}




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



TransferCellAnno <- function(Sobj_ref = NULL,
                         Sobj_query = NULL,
                         reference_cell_type_col = NULL,
                         reference_PCA_reduction = NULL,
                         reference_UMAP_reduction = NULL,
                         plot_results = NULL) {
dims <- c(1:20)  
DefaultAssay(Sobj_ref) <- "RNA"
DefaultAssay(Sobj_query) <- "RNA"
# Filter Sobj_ref to Sobj_query genes (5k)
Sobj_ref <- subset(Sobj_ref, features = rownames(Sobj_query))
# Run Seurat workflow 
Sobj_query <- NormalizeData(Sobj_query)
Sobj_query <- FindVariableFeatures(Sobj_query)
Sobj_query <- ScaleData(Sobj_query)
Sobj_query <- RunPCA(Sobj_query)
Sobj_query <- FindNeighbors(Sobj_query, dims = dims)
Sobj_query <- FindClusters(Sobj_query)

Sobj_ref <- NormalizeData(Sobj_ref)
Sobj_ref <- FindVariableFeatures(Sobj_ref)
Sobj_ref <- ScaleData(Sobj_ref)
Sobj_ref <- RunPCA(Sobj_ref)
Sobj_ref <- FindNeighbors(Sobj_ref, dims = dims)
Sobj_ref <- FindClusters(Sobj_ref)

# Projecting cell anno. to a query dataset
anchors <- FindTransferAnchors(reference = Sobj_ref, query = Sobj_query, dims = dims, reference.reduction = "pca", k.anchor = 5)

Sobj_ref <- RunUMAP(Sobj_ref, dims = dims, reduction = "pca", return.model = TRUE)

Sobj_query <- MapQuery(anchorset = anchors, reference = Sobj_ref, query = Sobj_query,
                       refdata = as.vector(Sobj_ref@meta.data[[reference_cell_type_col]]), 
                       reference.reduction = reference_PCA_reduction, 
                       reduction.model = reference_UMAP_reduction)

# Optionally plot results for the Sobj_query data
if (plot_results) {
  p1 <- DimPlot(Sobj_ref, reduction = "umap", group.by = reference_cell_type_col, label = TRUE, label.size = 5, repel = TRUE) + ggtitle("Reference annotations") + NoLegend()
  
  p2 <- DimPlot(Sobj_query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
                label.size = 5, repel = TRUE) + ggtitle("Query transferred labels") + NoLegend()
  
  suppressWarnings({
    print(p1 + p2)
  })
  
  
}

Sobj_list <- list(Sobj_query = Sobj_query,
                  Sobj_ref = Sobj_ref,
                  dimplot_query = p2,
                  dimplot_ref = p1)
return(Sobj_list)
}



### Cell Type Projection --------------------------

## Brain - Cells --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref_1 <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                              cell_type_pattern = "Brain_Neurons", 
                              metadata = metadata)
Sobj_ref_2 <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                              cell_type_pattern = "Brain_Microglia", 
                              metadata = metadata)
Sobj_ref <- Reduce(function(x, y) merge(x, y), list(Sobj_ref_1, Sobj_ref_2))
Sobj_ref[["RNA"]] <- JoinLayers(Sobj_ref[["RNA"]])


# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                Organ = "brain",
                                Cell_Nuc = "Cells",
                                High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                             Sobj_query = Sobj_query,
                             reference_cell_type_col = "cell_ontology_class",
                             plot_results = TRUE)

Sobj_anno_brain_cell <- Sobj_anno


## Heart - Cells --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                              cell_type_pattern = "Heart", 
                              metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "heart",
                                  Cell_Nuc = "Cells",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_heart_cell <- Sobj_anno


## Liver - Cells --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Liver", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "liver",
                                  Cell_Nuc = "Cells",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_liver_cell <- Sobj_anno




## Kidney - Cells --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Kidney", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "kidney",
                                  Cell_Nuc = "Cells",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_kidney_cell <- Sobj_anno




## Lung - Cells --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Lung", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "lung",
                                  Cell_Nuc = "Cells",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_lung_cell <- Sobj_anno







## Brain - Nucs --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref_1 <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                              cell_type_pattern = "Brain_Neurons", 
                              metadata = metadata)
Sobj_ref_2 <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                              cell_type_pattern = "Brain_Microglia", 
                              metadata = metadata)
Sobj_ref <- Reduce(function(x, y) merge(x, y), list(Sobj_ref_1, Sobj_ref_2))
Sobj_ref[["RNA"]] <- JoinLayers(Sobj_ref[["RNA"]])


# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "brain",
                                  Cell_Nuc = "Nucs",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_brain_nuc <- Sobj_anno


## Heart - Nucs --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Heart", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "heart",
                                  Cell_Nuc = "Nucs",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_heart_nuc <- Sobj_anno


## Liver - Nucs --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Liver", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "liver",
                                  Cell_Nuc = "Nucs",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_liver_nuc <- Sobj_anno




## Kidney - Nucs --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Kidney", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "kidney",
                                  Cell_Nuc = "Nucs",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_kidney_nuc <- Sobj_anno




## Lung - Nucs --------
# # var inputs 
# metadata 
# ref_data_dir 
# cell_type_pattern

Sobj_ref <- Read_Ref_SCdata(data_dir = ref_data_dir, 
                            cell_type_pattern = "Lung", 
                            metadata = metadata)

# # var inputs 
# Sobj_merg_CellNuc 
# Organ 
# Cell_or_Nuc
# High_Low_QC
Sobj_query <- Filter_Query_SCdata(Sobj = Sobj_merg_CellNuc,
                                  Organ = "lung",
                                  Cell_Nuc = "Nucs",
                                  High_Low_QC = "High quality cells")


# # var inputs 
# Sobj_ref 
# Sobj_query 
# reference_cell_type_col
# plot_results
Sobj_anno <- TransferCellAnno(Sobj_ref = Sobj_ref,
                              Sobj_query = Sobj_query,
                              reference_cell_type_col = "cell_ontology_class",
                              plot_results = TRUE)

Sobj_anno_lung_nuc <- Sobj_anno


# Define the lists for each organ's cell and nucleus
organ_cell_info <- list(
  Brain = Sobj_anno_brain_cell,
  Heart = Sobj_anno_heart_cell,
  Liver = Sobj_anno_liver_cell,
  Kidney = Sobj_anno_kidney_cell,
  Lung = Sobj_anno_lung_cell
)

organ_nuc_info <- list(
  Brain = Sobj_anno_brain_nuc,
  Heart = Sobj_anno_heart_nuc,
  Liver = Sobj_anno_liver_nuc,
  Kidney = Sobj_anno_kidney_nuc,
  Lung = Sobj_anno_lung_nuc
)

# Combine both lists into one
Sobj_organs_annotted <- list(Cells = organ_cell_info, Nuclei = organ_nuc_info)

# Save the combined information to an R data file
# saveRDS(Sobj_organs_annotted, file = "/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")


Sobj_organs_annotted_NOmix <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_organs_annotted_NOmix.rds")






# Define vectors for organs and modalities
Organs <- c("Heart", "Lung", "Liver", "Kidney", "Brain")  # Add more organs as needed
Modalities <- c("Cells", "Nuclei")    # Add more modalities as needed


for (Modality in Modalities) {
  for (Organ in Organs) {
    
    # Construct the filename based on the current Organ and Modality
    filename <- paste0("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/High_res/UMAPs_Mix_Highlight/", Organ, "_", Modality, "_Celltype.pdf")
    
    # Open a PNG device
    pdf(filename, width = 8, height = 6)
    
    # Create the UMAP plot with error handling
      DimPlot(Sobj_organs_annotted_NOmix[[Modality]][[Organ]]$Sobj_query,
              reduction = "ref.umap", 
              group.by = "predicted.id", 
              label = FALSE, 
              label.size = 4, 
              pt.size = 1,
              repel = FALSE) +
        #ggtitle(paste0(Organ, "_", Modality)) +
        ggtitle("") +
        theme_void() +  
        theme(panel.border = element_blank(),  
              plot.background = element_blank(),
              legend.text = element_text(size = 15),  
              legend.title = element_text(size = 25))  
      
      # Close the PNG device
      dev.off()
  }
}



## Mix cells 
# Ran as /mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/scripts/Mix_cells_celltype_NOsketch_YESjoinlayers.R
Sobj_Cells_Mix_annotted <- readRDS("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/objects/Sobj_Cells_Mix_annotted_NOorgans.rds")


# filename <- paste0("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/High_res/UMAPs_Mix_Highlight/", "Mix", "_", "Cells", "_Celltype.pdf")

# Open a PNG device
pdf(filename, width = 8, height = 6)

# Create the UMAP plot with error handling
DimPlot(Sobj_Cells_Mix_annotted,    ## no good - too much label for a UMAP
        reduction = "ref.umap", 
        group.by = "predicted.id",
        label = FALSE, 
        label.size = 4, 
        pt.size = 1,
        repel = FALSE) +
  #ggtitle(paste0(Organ, "_", Modality)) +
  ggtitle("") +
  theme_void() +  
  theme(panel.border = element_blank(),  
        plot.background = element_blank(),
        legend.text = element_text(size = 15),  
        legend.title = element_text(size = 25))  

# Close the PNG device
dev.off()



# Not good - too many labels
filename <- paste0("/mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/assets/High_res/UMAPs_Mix_Highlight/", "Mix", "_", "Cells", "_Celltype_Barplot.pdf")

# Open a PNG device
pdf(filename, width = 8, height = 6)
library(dittoSeq)
dittoBarPlot(Sobj_Cells_Mix_annotted, 
             "predicted.id", 
             group.by = "Tissue_type") +
  labs(title = "Cell Type vs Organ projection (Cells)",
       x = "Organ",
       y = "Percent of Cells") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
    # , legend.position = "none")
  )

dev.off()






## Mix nucs
# Ran as /mnt/scratch1/maycon/STAMP.git/Paper_QC/round_1/scripts/Mix_nucs_celltype_NOsketch_YESjoinlayers.R











