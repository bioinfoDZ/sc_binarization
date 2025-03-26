# Multigrate Integration

This directory contains scripts for multimodal single-cell data integration using the multigrate package. 

## Contents

- **sim_multigrate.py**: Core integration script that performs multimodal data integration and clustering.
  - Loads RNA and ATAC-seq data
  - Integrates data using MultiVAE model
  - Performs clustering (Leiden algorithm)
  - Evaluates clustering quality against ground truth
  - Generates visualizations and metrics

## Dependencies

- numpy
- scanpy
- multigrate
- muon
- pandas
- scikit-learn
- anndata
- matplotlib

