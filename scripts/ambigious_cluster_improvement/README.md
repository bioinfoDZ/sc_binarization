# Ambiguous Clusters Analysis

This directory contains scripts for analyzing and visualizing ambiguous cell clusters in single-cell data, particularly focusing on monocyte subtypes in PBMC datasets.

## Contents

- **01_process_ambiguous_clusters.py**: Core analysis script that performs binary parsing of single-cell data and evaluates clustering accuracy.
  - Calculates clustering metrics for different ratios of genes/peaks
  - Generates cell type and cluster statistics
  - Outputs processed data for downstream visualization

- **02_plot_atac_purity.py**: Visualizes the effect of varying gene/peak ratios on clustering purity.
  - Focuses on ATAC-based standards
  - Generates line plots showing purity across different percentages


- **03_plot_rna_atac_comparison.py**: Creates comparison visualizations for RNA and ATAC standards.
  - Analyzes ambiguous cell types (intermediate monocytes and CD14 monocytes)
  - Compares clustering performance using different annotation standards
  