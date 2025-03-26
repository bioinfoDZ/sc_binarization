# Simulation

This directory contains scripts for generating synthetic single-cell multi-modal data for benchmarking and testing integration methods.

## Contents

    simulation.R: R script for simulating paired RNA-seq and ATAC-seq data.
        Uses the scMultiSim package to generate realistic synthetic data
        Creates data with known ground truth populations
        Simulates gene regulatory networks (GRNs) with controlled parameters
        Outputs data in 10X Genomics compatible format

## Dependencies

    R (>= 4.0.0)
    Required R packages:
        scMultiSim
        Matrix
        optparse
        ape

## Usage

Run the simulation script with custom parameters:


Rscript simulation.R --index 1 --num_cells 10000 --num_genes 2500 --output_dir sim_results

## Parameters

    --index: Integer identifier for the simulation run (default: 1)
    --num_cells: Number of cells to simulate (default: 10000)
    --num_genes: Number of genes to simulate (default: 2500)
    --cif_sigma: Cell-intrinsic factor sigma value (default: 0.5)
    --diff_cif_fraction: Differentiation cell-intrinsic factor fraction (default: 0.2)
    --intrinsic_noise: Intrinsic noise level (default: 0.2)
    --output_dir: Directory to save simulation outputs

Output Files

The script generates several output files in the specified output directory:

    matrix.mtx: Sparse matrix in Matrix Market format containing combined RNA and ATAC-seq counts
    features.tsv: Tab-separated file containing feature information (genes and regions)
    barcodes.tsv: Tab-separated file containing cell barcodes with population annotations
    [index]_index_tree.png: Visualization of the phylogenetic tree used for simulation
    [index]_newick_tree.txt: Newick format representation of the tree structure
