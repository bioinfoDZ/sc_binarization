# Single Cell Binarized Integration

This repository contains tools and workflows for single-cell data binarization and integration across multiple modalities, as described in our manuscript on single cell binary integration methods.

## Project Description

Our approach transforms continuous gene expression data into binary states, facilitating novel integration methods for multi-modal single-cell datasets. This repository provides workflows for data binarization, integration, and downstream analysis.

## Repository Structure

    data/: Example datasets used in the manuscript
    scripts/:
        Simulation scripts for generating synthetic single-cell data
        Integration software implementation scripts
    jupyter_notebooks/: Jupyter notebooks demonstrating complete analysis workflows
    r_scripts/: R scripts for dataset-specific analyses
        Analysis of multiomics cortex data (Zhu et al., Sci Adv 2023)
        Bone marrow stroma scRNA-seq analysis (Baryawno et al., Cell 2019)
        Feature inclusion sweep visualizations for murine breast dataset

## Requirements

Python Environment

    Python 3.11.4
    Scanpy 1.9.3
    Additional dependencies listed in environment.yml

R Environment

    R 4.1.0
    Seurat 5.1.0
    Signac 1.11.0
