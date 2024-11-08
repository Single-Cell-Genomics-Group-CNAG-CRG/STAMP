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


# Load necessary package
library(msigdbr)

# Define the pathway names
pathways <- c("HALLMARK_APOPTOSIS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_HYPOXIA", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_P53_PATHWAY", "BIOCARTA_STEM_PATHWAY",
              "BIOCARTA_LYMPHOCYTE_PATHWAY")

# Initialize an empty list to store the vectors
pathway_genes <- list()
# Loop through each pathway and store the genes in a named vector within the list
for (pathway in pathways) {
  # Filter genes for the specific pathway
  genes <- msigdbr(species = "Homo sapiens", category = ifelse(grepl("BIOCARTA", pathway), "C2", "H")) %>%
    dplyr::filter(gs_name == pathway) %>%
    dplyr::pull(gene_symbol)
  
  # Store genes in the list with the pathway name as the list element name
  pathway_genes[[pathway]] <- genes
}
# Example: Access genes in HALLMARK_APOPTOSIS
pathway_genes$HALLMARK_APOPTOSIS

library(AUCell)
scoring <- AUCell_run(counts(sce), 
                        geneSets = pathway_genes,
                        BPPARAM = bp)


