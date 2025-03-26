# Muon Integration

This directory contains scripts for multimodal single-cell data integration using the Muon package, specifically implementing MultiOmics Factor Analysis (MOFA).

## Contents

- **sim_muon.py**: Integration script that performs multimodal data integration using MOFA.
  - Loads RNA and ATAC-seq data from separate h5ad files
  - Initializes a pyMOFA model through the omicverse package
  - Preprocesses the data for integration
  - Runs factor analysis with GPU acceleration
  - Saves the integrated model to an HDF5 file

## Dependencies

- muon
- scanpy
- cupy
- omicverse
