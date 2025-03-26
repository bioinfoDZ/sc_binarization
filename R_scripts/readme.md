# Single-Cell Data Binarization and Analysis

This folder contains R scripts for the analysis of single-cell sequencing data, including RNA binarization, clustering, and visualization across multiple datasets.


## Scripts Overview

10x Multiomics Cortex Dataset (Zhu et al., Sci Adv 2023)

    cortex_script.R: Analysis pipeline for 10x multiomics RNA and ATAC-seq data, including RNA count binarization and clustering
    cortex_figures.R: Generation of manuscript figures for the reanalysis of the cortex dataset

Bone Marrow Stroma Dataset (Baryawno et al., Cell 2019)

    readin_bmstroma.R: Data parsing and preprocessing of bone marrow stroma (bone cell) scRNA-seq data
    final_boneatlas_script.R: Complete analysis workflow including binarization and clustering of bone marrow stroma cells
    final_boneatlas_figures.R: Production of manuscript figures for the bone atlas reanalysis

Additional Analysis

    murinebreast_hvfsweep_figures.R: Preparation of violin plots summarizing feature inclusion sweep for murine breast dataset (Figure 4B)

Requirements

    R (version 4.0 or higher recommended)
    Required R packages:
        Seurat
        ggplot2
        dplyr
        Matrix
        patchwork
        umap
        additional packages as imported by individual scripts
