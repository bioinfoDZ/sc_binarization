#!/usr/bin/env Rscript

library(scMultiSim)
library(Matrix)
library(optparse)
library(ape)


# Setup command-line argument parsing
option_list <- list(
  make_option(c("--index"), type = "integer", default = 1, help = "Index", metavar = "integer"),
  make_option(c("--num_cells"), type = "integer", default = 10000, help = "Number of cells", metavar = "integer"),
  make_option(c("--num_genes"), type = "integer", default = 2500, help = "Number of genes", metavar = "integer"),
  make_option(c("--cif_sigma"), type = "numeric", default = 0.5, help = "CIF sigma value", metavar = "number"),
  make_option(c("--diff_cif_fraction"), type = "numeric", default = 0.2, help = "Diff CIF fraction", metavar = "number"),
  make_option(c("--intrinsic_noise"), type = "numeric", default = 0.2, help = "Intrinsic noise level", metavar = "number"),
  make_option(c("--output_dir"), type = "character", default = "sim_results_new", help = "Output directory", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load the provided GRN parameters
data(GRN_params_1139)
GRN_params <- GRN_params_1139

rand_tree = Phyla5() #ape::rtree(7)
png(paste0(opt$index,"_index_tree.png"))
plot(rand_tree, show.tip.label = TRUE, )
# Add edge lengths (distances) to the plot
#edgelabels(round(rand_tree$edge.length, 2), bg="white")
dev.off()
write.tree(rand_tree, file = paste0(opt$index,"_newick_tree.txt"), append = FALSE,
           digits = 10, tree.names = FALSE)


arith_progression_sum <- function(total_sum, start_value, num_values) {
  common_difference <- (2 * (total_sum - start_value * num_values)) / (num_values * (num_values - 1))
  sequence <- start_value + (0:(num_values - 1)) * common_difference
  integer_sequence <- round(sequence)
  return(integer_sequence)
}


total_sum <- 5000

start_value <- 300

num_values <- 5

cluster_sizes <- arith_progression_sum(total_sum, start_value, num_values)
print(cluster_sizes)
print(class(cluster_sizes))
print(sum(cluster_sizes))



# Define simulation options
options <- list(
  rand.seed = 0,
  GRN = GRN_params,
  num.cells = opt$num_cells,
  num.genes = opt$num_genes,
  cif.sigma = opt$cif_sigma,
  tree = rand_tree,
  diff.cif.fraction = opt$diff_cif_fraction,
  do.velocity = FALSE,
  discrete.cif = TRUE,
  intrinsic.noise = opt$intrinsic_noise,
  discrete.min.pop.size = cluster_sizes[1],
  discrete.min.pop.index = 1,
  discrete.pop.size = as.integer(cluster_sizes) # c(250, 500, 750, 1000, 1250, 1500, 1750)
)

print('running sim')
# Run simulation
res <- sim_true_counts(options)
print('done')


# Extract and process results
cnt <- res$counts
atac_cnt <- res$atacseq_data
mt <- res$cell_meta
r2g <- res$region_to_gene

# Generate feature identifiers and cell names
genes <- paste0("gene", seq_len(nrow(cnt)))
regions <- paste0("region", seq_len(nrow(atac_cnt)))
cells <- paste0("cell", seq_len(ncol(cnt)), "_type", mt$pop)

# Ensure the specific output directory exists
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Concatenate matrices vertically
combined_matrix <- rbind(Matrix(cnt, sparse = TRUE), Matrix(atac_cnt, sparse = TRUE))

# Combined feature identifiers and types
combined_features <- c(genes, regions)
feature_types <- c(rep("Gene Expression", length(genes)), rep("ATAC", length(regions)))

# Write the concatenated matrix to Matrix Market format
writeMM(combined_matrix, file.path(opt$output_dir, "matrix.mtx"))

# Write combined features to TSV, including feature type
write.table(cbind(combined_features, combined_features, feature_types), file.path(opt$output_dir, "features.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Write cell barcodes to TSV
write.table(cells, file.path(opt$output_dir, "barcodes.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Compress the files
system(paste("gzip -f", file.path(opt$output_dir, "matrix.mtx")))
system(paste("gzip -f", file.path(opt$output_dir, "features.tsv")))
system(paste("gzip -f", file.path(opt$output_dir, "barcodes.tsv")))
