# STAMP: Single-Cell Transcriptomics Analysis and Multimodal Profiling Through Imaging

STAMP (Single-Cell Transcriptomics Analysis and Multimodal Profiling) is a scalable, cost-effective approach for single-cell profiling that eliminates sequencing costs by integrating spatial transcriptomics and proteomics imaging. By immobilizing cells onto slides, STAMP enables high-throughput single- and multimodal (RNA, protein, and H&E) profiling while preserving cellular structure. Its flexible format supports large-scale studies across PBMCs, cell lines, stem cells, dissociated tissues, and FFPE samples. STAMP facilitates high-throughput immunophenotyping, rare cell detection, and perturbation studies, demonstrated across over 10 million cells and 6 billion transcripts. This strategy revolutionizes cellular profiling, making it more accessible and scalable for research and clinical applications.

This repository contains all the scripts, notebooks and reports to reproduce all analysis from our manuscript entitled STAMP: Single-Cell Transcriptomics and Multimodal Profiling Through Imaging.


## Technologies Used
- **Spatial Transcriptomics Platforms**:
  - NanoString CosMx
  - 10X Genomics Xenium
- **Single Cell RNA sequencing Platforms**:
  - 10X Genomics Flex
  - 10X Genomics 5'
  - 10X Genomics 3'
 
## Repository Structure

- **misc/**
  - `BIN.R`: Contains functions used throughout the analyses.
  - `paths`: Contains various paths used in the project.
- **stamp_1/**, **stamp_2/**, ..., **stamp_n/**:
  - Each folder corresponds to a sample and contains the code specific to that sample.
  - **Inside each sample folder:**
    - **Preparation/**
      - Scripts for loading the expression matrix and creating `SingleCellExperiment` objects.
    - **QC/**
      - Scripts used to perform quality control.
    - **Analysis/**
      - `PreProc.R`: Performs preprocessing steps like normalization, feature selection, and dimensionality reduction.
      - `Clust.R`: Performs clustering analysis.
    - **Lvl1/** and **Lvl2/**:
      - Scripts used for different rounds of cell population annotation.
    - **Sub-sample Folders**:
      - For samples containing multiple sub-samples (e.g., multiplexed slides), there are folders like **Clines/** and **PBMC/** (e.g., in `stamp_3/`), which contain analyses of the respective sub-samples.

## R Version

Most of this code was run using **R version 4.4.1** in a MacBook Pro M3. Computational heavy tasks, such as alignment of the scRNAseq FASTQs to the reference genome, were performed in a HPC cluster.
---

*For more details, please refer to the project's paper or contact the authors.*

