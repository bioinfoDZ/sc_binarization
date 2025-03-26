# Single-Cell Bayesian Factor Analysis (scBFA)

This directory contains R scripts for performing Bayesian Factor Analysis on single-cell data. The implementation allows for dimension reduction and feature extraction from high-dimensional single-cell datasets.

## Contents

    build_model.R: Script for constructing and training the Bayesian Factor Analysis model.
        Sets up model parameters
        Fits BFA model to single-cell expression data
        Extracts latent factors and loadings
        Saves model output for downstream analysis

    analyze_models.R: Script for evaluating and visualizing BFA model results.
        Loads trained models
        Generates diagnostic plots
        Computes quality metrics
        Compares different model parameterizations
        Creates visualizations of latent factors
