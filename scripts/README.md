# Single-Cell Data Binarization and Integration

This folder contains scripts for single-cell data binarization, integration, and analysis across multiple modalities. The project implements several methods for multimodal data integration and factor analysis.

## Folder Structure

    simulation/: Scripts for generating synthetic single-cell multi-modal data
        Realistic simulations of RNA-seq and ATAC-seq data with ground truth populations
        Customizable parameters for various simulation scenarios

    scbfa/: Single-Cell Bayesian Factor Analysis implementation
        Model building and training scripts
        Analysis and visualization tools for BFA results

    muon/: Muon-based integration of multi-modal data
        Implementation of MOFA (Multi-Omics Factor Analysis)
        Processing and analysis of integrated data

    momat/: scMOMAT integration methods
        Maximum-margin optimization for multi-omics alignment
        Clustering and visualization of integrated results

    multigrate/: multigrate integration 
        Constructing weighted nearest-neighbor graphs for each modality
        Computing joint feature spaces that preserve modality-specific signals

    ambigious_cluster_improvement/: analyzing ambiguous cell clusters
        Identifies cells with uncertain classification using ensemble methods
        Quantifies the degree of ambiguity for each cell