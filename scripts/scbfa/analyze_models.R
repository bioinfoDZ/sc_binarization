#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(scBFA)
    library(umap)
    library(ggplot2)
    library(mclust)  # For ARI
    library(aricode) # For AMI, NMI
})

# Load sampled cell barcodes and cell types
sampled_barcodes <- read.csv("outputs/sampled_cell_barcodes.csv", stringsAsFactors = FALSE)[[1]]
annot_data <- read.csv("pbmc10k_celltypes.csv", row.names = 1)
sampled_types <- annot_data[sampled_barcodes, "atac.celltype"]

# Filter to only cells with valid cell type annotations and log the count
valid_indices <- !is.na(sampled_types)
sampled_types <- factor(sampled_types[valid_indices])  # Ensure no NA in cell types
sampled_barcodes <- sampled_barcodes[valid_indices]
num_valid_cells <- length(sampled_types)
cat(paste("Total valid cells used in analysis:", num_valid_cells, "\n"))

# Function to analyze a model, plot UMAP, factor plots, and calculate clustering metrics
analyze_model <- function(model_file, sampled_types, plot_prefix) {
    # Load the model
    model <- readRDS(model_file)
    
    # Extract model factors (cell embeddings in low-dimensional space)
    model_factors <- model$ZZ
    model_factors <- model_factors[valid_indices, , drop = FALSE]  # Align with valid cell barcodes
    
    # Run UMAP on the factors
    umap_result <- umap::umap(model_factors)
    
    # Perform K-means clustering
    kmeans_result <- kmeans(model_factors, centers = length(unique(sampled_types)))
    
    # Calculate clustering metrics, handling potential NA issues
    metrics <- list(
        ARI = mclust::adjustedRandIndex(sampled_types, kmeans_result$cluster),
        AMI = if (!any(is.na(sampled_types))) aricode::AMI(sampled_types, kmeans_result$cluster) else NA,
        MI = if (!any(is.na(sampled_types))) aricode::NMI(sampled_types, kmeans_result$cluster) else NA
    )
    
    # Prepare data frame for plotting
    plot_df <- data.frame(
        UMAP1 = umap_result$layout[, 1],
        UMAP2 = umap_result$layout[, 2],
        Factor1 = model_factors[, 1],  # Factor 1
        Factor2 = model_factors[, 2],  # Factor 2
        CellType = factor(sampled_types),
        Cluster = factor(kmeans_result$cluster)
    )
    
    # Calculate axis limits for tighter plots
    umap_x_limits <- range(plot_df$UMAP1) + c(-0.1, 0.1) * diff(range(plot_df$UMAP1))
    umap_y_limits <- range(plot_df$UMAP2) + c(-0.1, 0.1) * diff(range(plot_df$UMAP2))
    factor_x_limits <- range(plot_df$Factor1) + c(-0.1, 0.1) * diff(range(plot_df$Factor1))
    factor_y_limits <- range(plot_df$Factor2) + c(-0.1, 0.1) * diff(range(plot_df$Factor2))

    # File paths for outputs
    file_paths <- list(
        celltype_umap = paste0("outputs/plots/", plot_prefix, "_umap_celltype_plot.png"),
        cluster_umap = paste0("outputs/plots/", plot_prefix, "_umap_cluster_plot.png"),
        celltype_factor = paste0("outputs/plots/", plot_prefix, "_factor12_celltype_plot.png"),
        cluster_factor = paste0("outputs/plots/", plot_prefix, "_factor12_cluster_plot.png"),
        metrics = paste0("outputs/metrics/", plot_prefix, "_clustering_metrics.txt")
    )
    
    # Delete existing files if they exist
    lapply(file_paths, function(fp) if (file.exists(fp)) file.remove(fp))
    
    # Generate UMAP plot with cell types
    png(filename = file_paths$celltype_umap, width = 1000, height = 1000, res=120)
    p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
        geom_point(size = 2, alpha = 1) +
        # Remove axis ticks and text
        scale_x_continuous(labels = NULL, breaks = NULL) +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        ggtitle(paste("UMAP of", plot_prefix, "Factors by Cell Type")) +
        xlim(umap_x_limits) + ylim(umap_y_limits) +
        theme_minimal(base_size=16) +
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16)
        )
    print(p1)
    dev.off()
    
    # Generate UMAP plot with clusters and annotated metrics
    png(filename = file_paths$cluster_umap, width = 1000, height = 1000, res=120)
    p2 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 2, alpha = 1) +
                # Remove axis ticks and text
        scale_x_continuous(labels = NULL, breaks = NULL) +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        ggtitle(paste("UMAP of", plot_prefix, "Factors with Clusters")) +
        xlim(umap_x_limits) + ylim(umap_y_limits) +
        theme_minimal(base_size=16) +
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16)
        )
    print(p2)
    dev.off()
    
    # Generate Factor 1 vs Factor 2 plot with cell types
    png(filename = file_paths$celltype_factor, width = 1000, height = 1000, res=120)
    p3 <- ggplot(plot_df, aes(x = Factor1, y = Factor2, color = CellType)) +
        geom_point(size = 2, alpha = 1) +
        # Change axis labels to tSNE
        labs(x = "tSNE axis 1", y = "tSNE axis 2") +
        # Remove axis ticks and text
        scale_x_continuous(labels = NULL, breaks = NULL) +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        ggtitle(paste("tSNE of", plot_prefix, "by Cell Type")) +
        xlim(factor_x_limits) + ylim(factor_y_limits) +
        theme_minimal(base_size=16) +
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        # Explicitly remove axis text
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16)
        )
    print(p3)
    dev.off()
    
    # Generate Factor 1 vs Factor 2 plot with clusters
    png(filename = file_paths$cluster_factor, width = 1000, height = 1000, res=120)
    p4 <- ggplot(plot_df, aes(x = Factor1, y = Factor2, color = Cluster)) +
        geom_point(size = 2, alpha = 1) +
        # Change axis labels to tSNE
        labs(x = "tSNE axis 1", y = "tSNE axis 2") +
        # Remove axis ticks and text
        scale_x_continuous(labels = NULL, breaks = NULL) +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        ggtitle(paste("tSNE of", plot_prefix, "with Clusters")) +
        xlim(factor_x_limits) + ylim(factor_y_limits) +
        theme_minimal(base_size=16) +
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16)
        )
    print(p4)
    dev.off()
    
    # Save clustering metrics to a text file
    write.table(data.frame(
        Metric = names(metrics),
        Value = unlist(metrics)
    ), 
    file_paths$metrics,
    row.names = FALSE, quote = FALSE, sep = "\t")
    
    print(paste(plot_prefix, "model analysis complete! Total valid cells used:", num_valid_cells))
}

# Run analysis on each model
analyze_model("outputs/full_model.rds", sampled_types, "full_rna_atac")
analyze_model("outputs/gene_only_model.rds", sampled_types, "gene_only")
analyze_model("outputs/peak_only_model.rds", sampled_types, "peak_only")

print("All analyses complete.")