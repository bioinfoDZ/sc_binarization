#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
model_type <- args[1]  # "full", "gene_only", or "peak_only"
output_model <- args[2]  # Path to save the model

suppressPackageStartupMessages({
    library(Seurat)
    library(scBFA)
})

# Load RNA+ATAC data and sampled cell barcodes
data_list <- Read10X_h5("10k_pbmcs/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_barcodes <- colnames(data_list[["Gene Expression"]])
atac_barcodes <- colnames(data_list[["Peaks"]])
sampled_barcodes <- read.csv("outputs/sampled_cell_barcodes.csv", stringsAsFactors = FALSE)[[1]]

# Check that all sampled_barcodes are valid
valid_barcodes <- intersect(rna_barcodes, atac_barcodes)  # Only barcodes present in both RNA and ATAC
missing_barcodes <- setdiff(sampled_barcodes, valid_barcodes)
if (length(missing_barcodes) > 0) {
    stop("The following barcodes are missing from the data: ", paste(missing_barcodes, collapse = ", "))
}

# Prepare data based on model type
if (model_type == "full") {
    # Combine RNA and ATAC matrices for the full model
    rna_matrix <- as.matrix(data_list[["Gene Expression"]][, sampled_barcodes, drop = FALSE])
    atac_matrix <- as.matrix(data_list[["Peaks"]][, sampled_barcodes, drop = FALSE])
    combined_matrix <- rbind(rna_matrix, atac_matrix)
    colnames(combined_matrix) <- paste0("GeneExpr_", colnames(combined_matrix))  # Label all as "Gene Expression"
    scData <- combined_matrix
} else if (model_type == "gene_only") {
    # Use only RNA data and label as "Gene Expression"
    rna_matrix <- as.matrix(data_list[["Gene Expression"]][, sampled_barcodes, drop = FALSE])
    colnames(rna_matrix) <- paste0("GeneExpr_", colnames(rna_matrix))
    scData <- rna_matrix
} else if (model_type == "peak_only") {
    # Use only ATAC data and label as "Gene Expression"
    atac_matrix <- as.matrix(data_list[["Peaks"]][, sampled_barcodes, drop = FALSE])
    colnames(atac_matrix) <- paste0("GeneExpr_", colnames(atac_matrix))
    scData <- atac_matrix
} else {
    stop("Unknown model type: ", model_type)
}

# Build model with 6 factors
model <- scBFA(
    scData = scData,  # Ensure input is a single matrix
    numFactors = 6,
    method = "CG",
    maxit = 300
)

# Save model
saveRDS(model, output_model)
