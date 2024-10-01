# STAMP: Single-Cell Transcriptomics Analysis and Multimodal Profiling through Imaging

This repository contains the code used in the **STAMP** project, focusing on spatial transcriptomics analysis and multiplexing profiling of single cells. The code was utilized to analyze the data presented in the corresponding paper.

## Technologies Used
- **Spatial Transcriptomics Platforms**:
  - NanoString CosMx
  - 10X Genomics Xenium
- **Single Cell RNA sequencing Platforms**:
  - 10X Genomics Flex
 
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

This code was run using **R version 4.4.1**.
---

*For more details, please refer to the project's paper or contact the authors.*

