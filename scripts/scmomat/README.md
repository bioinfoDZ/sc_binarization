# scMOMAT Integration

This directory contains scripts for multimodal single-cell data integration using the scMOMAT (Single-Cell Multi-Omics alignment with Maximum-margin Optimization) method.

## Contents

    sim_momat.py: Integration script that performs multimodal data integration and analysis.
        Loads RNA and ATAC-seq data
        Preprocesses single-cell omics data for integration
        Trains a scMOMAT model with factor extraction
        Performs dimensionality reduction (PCA and UMAP)
        Conducts clustering using the Leiden algorithm
        Calculates quality metrics for clusters
        Generates visualizations of integrated data

Dependencies

    scmomat
    scanpy
    muon
    numpy
    pandas
    scikit-learn
    matplotlib
    umap-learn
    scipy
    torch
