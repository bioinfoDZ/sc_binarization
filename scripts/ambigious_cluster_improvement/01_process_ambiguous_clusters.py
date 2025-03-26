import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score
from muon import atac as ac
import scipy.sparse
import os
from datetime import datetime
import sys
import warnings
import h5py
warnings.filterwarnings('ignore')

# File structure setup
BASE_DIR = os.getcwd()
INPUT_DIR = os.path.join(BASE_DIR, "input_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
PLOT_DIR = os.path.join(OUTPUT_DIR, "plots")
DATA_DIR = os.path.join(OUTPUT_DIR, "processed_data")
RESULTS_DIR = os.path.join(OUTPUT_DIR, "results")
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")

# Create all necessary directories
for d in [INPUT_DIR, OUTPUT_DIR, PLOT_DIR, DATA_DIR, RESULTS_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

def setup_logging():
    """Setup logging to both file and console"""
    log_file = os.path.join(LOG_DIR, f"analysis_{datetime.now():%Y%m%d_%H%M%S}.log")
    
    class Logger(object):
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w")
        
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
            self.flush()
            
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    
    sys.stdout = Logger(log_file)

def check_file_integrity(filename):
    """Check if h5ad file is readable"""
    try:
        with h5py.File(filename, 'r') as f:
            required_groups = ['X', 'obs', 'var']
            for group in required_groups:
                if group not in f:
                    return False, f"Missing required group: {group}"
        return True, "File is valid"
    except Exception as e:
        return False, str(e)

def safe_load_h5ad(filename):
    """Safely load h5ad file with error handling"""
    print(f"Attempting to load: {filename}")
    is_valid, message = check_file_integrity(filename)
    
    if not is_valid:
        raise ValueError(f"File integrity check failed: {message}")
    
    try:
        return sc.read_h5ad(filename)
    except Exception as e:
        raise ValueError(f"Error loading h5ad file: {str(e)}")

def calculate_clustering_accuracy(adata, cell_type_key='rna:celltype', cluster_key='Binary_LSI_cluster'):
    """Calculate detailed clustering metrics"""
    cluster_stats = {}
    cell_type_stats = {}
    
    # Calculate stats for each cluster
    for cluster in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_cells = adata.obs[cluster_mask]
        
        type_counts = cluster_cells[cell_type_key].value_counts()
        dominant_type = type_counts.index[0]
        cluster_size = len(cluster_cells)
        purity = type_counts[0] / cluster_size
        
        cluster_stats[cluster] = {
            'dominant_type': dominant_type,
            'purity': purity,
            'size': cluster_size,
            'composition': type_counts.to_dict()
        }
    
    # Calculate stats for each cell type
    for cell_type in adata.obs[cell_type_key].unique():
        type_clusters = [
            cluster for cluster, stats in cluster_stats.items()
            if stats['dominant_type'] == cell_type
        ]
        
        if type_clusters:
            total_cells = sum(cluster_stats[c]['size'] for c in type_clusters)
            weighted_purity = sum(
                cluster_stats[c]['purity'] * cluster_stats[c]['size']
                for c in type_clusters
            ) / total_cells
            fragmentation = len(type_clusters)
        else:
            weighted_purity = 0
            fragmentation = 0
        
        cell_type_stats[cell_type] = {
            'weighted_purity': weighted_purity,
            'n_clusters': fragmentation,
            'total_cells': sum(adata.obs[cell_type_key] == cell_type)
        }
    
    return cell_type_stats, cluster_stats

def plot_clustering_metrics(results_df, iteration, output_dir):
    """Create visualizations for clustering metrics"""
    # Purity by cell type
    plt.figure(figsize=(12, 6))
    sns.barplot(data=results_df, x=results_df.index, y='weighted_purity')
    plt.xticks(rotation=45, ha='right')
    plt.title(f'Cluster Purity by Cell Type (Peak Ratio: {iteration:.1f})')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'purity_ratio_{iteration:.1f}.png'))
    plt.close()
    
    # Fragmentation by cell type
    plt.figure(figsize=(12, 6))
    sns.barplot(data=results_df, x=results_df.index, y='n_clusters')
    plt.xticks(rotation=45, ha='right')
    plt.title(f'Cluster Fragmentation by Cell Type (Peak Ratio: {iteration:.1f})')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'fragmentation_ratio_{iteration:.1f}.png'))
    plt.close()

def parser_bin_orig(adata, resolution=1, n_top_genes=25000):
    """Binary parser function"""
    adata_bin = adata.copy()
    
    # Convert X to proper array format before TFIDF
    if scipy.sparse.issparse(adata_bin.X):
        adata_bin.X = np.asarray(adata_bin.X.todense())
    else:
        adata_bin.X = np.asarray(adata_bin.X)
        
    ac.pp.tfidf(adata_bin, scale_factor=1e4)
    adata_bin.raw = adata_bin
    
    sc.experimental.pp.recipe_pearson_residuals(
        adata_bin, 
        n_top_genes=n_top_genes, 
        inplace=True
    )
    
    ac.tl.lsi(adata_bin)
    lsi = np.asarray(adata_bin.obsm['X_lsi'])
    varm = np.asarray(adata_bin.varm["LSI"])
    stdev = np.asarray(adata_bin.uns["lsi"]["stdev"])
    
    adata_bin.obsm['X_lsi'] = lsi[:, 1:]
    adata_bin.varm["LSI"] = varm[:, 1:]
    adata_bin.uns["lsi"]["stdev"] = stdev[1:]
    adata_bin.obsm['X_pca'] = adata_bin.obsm['X_lsi']
    
    sc.pp.neighbors(adata_bin, use_rep='X_lsi', n_pcs=15)
    sc.tl.umap(adata_bin)
    sc.tl.leiden(adata_bin, resolution=resolution, key_added='Binary_LSI_cluster')
    
    return adata_bin

def analyze_feature_ratio(adata, hvg_genes, hvg_peaks, peak_ratio=0.0):
    """Analyze clustering with different ratios of HVG genes and peaks"""
    adata_copy = adata.copy()
    
    n_peaks_to_add = int(len(hvg_peaks) * peak_ratio)
    selected_peaks = hvg_peaks[:n_peaks_to_add]
    
    print(f"Using all {len(hvg_genes)} HV genes + {n_peaks_to_add} peaks")
    
    adata_copy.var['highly_variable'] = False
    adata_copy.var.loc[hvg_genes, 'highly_variable'] = True
    adata_copy.var.loc[selected_peaks, 'highly_variable'] = True
    
    adata_processed = parser_bin_orig(adata_copy)
    
    # Remove large attributes before saving
    if 'pearson_residuals_df' in adata_processed.uns:
        del adata_processed.uns['pearson_residuals_df']
    
    cell_type_stats, cluster_stats = calculate_clustering_accuracy(adata_processed)
    
    return cell_type_stats, cluster_stats, adata_processed

def main():
    setup_logging()
    print(f"Starting analysis at {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Working directory: {BASE_DIR}")
    print(f"NumPy version: {np.__version__}")
    
    # Input file paths
    adata_file = os.path.join(INPUT_DIR, "pbmc10k_merge.h5ad")
    metadata_file = os.path.join(INPUT_DIR, "pbmc10k_celltypes.csv")
    
    # Check if input files exist
    if not os.path.exists(adata_file):
        raise FileNotFoundError(f"AnnData file not found: {adata_file}")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")
    
    # Load data with error handling
    try:
        print("Loading data...")
        adata = safe_load_h5ad(adata_file)
        metadata = pd.read_csv(metadata_file, index_col=0)
        adata.obs = adata.obs.join(metadata[['rna:celltype']])
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        sys.exit(1)
    
    # Initial processing
    print("Initial processing...")
    adata_bin = parser_bin_orig(adata)
    hvg_mask = adata_bin.var['highly_variable']
    
    # Get HV genes and peaks
    hvg_genes = adata.var.index[
        (hvg_mask) & (adata.var['feature_types'] == 'Gene Expression')
    ]
    hvg_peaks = adata.var.index[
        (hvg_mask) & (adata.var['feature_types'] == 'Peaks')
    ]
    
    print(f"Total HV features: {sum(hvg_mask)}")
    print(f"HV genes: {len(hvg_genes)}")
    print(f"HV peaks: {len(hvg_peaks)}")
    
    # Analysis for each ratio
    peak_ratios = np.arange(0, 1.1, 0.1)
    overall_results = []
    
    for ratio in peak_ratios:
        print(f"\nAnalyzing peak ratio: {ratio:.1f}")
        
        # Check for existing results
        iteration_file = os.path.join(DATA_DIR, f'iteration_{ratio:.1f}.h5ad')
        if os.path.exists(iteration_file):
            print(f"Loading existing results for ratio {ratio:.1f}")
            adata_processed = safe_load_h5ad(iteration_file)
            cell_type_stats, cluster_stats = calculate_clustering_accuracy(adata_processed)
        else:
            cell_type_stats, cluster_stats, adata_processed = analyze_feature_ratio(
                adata, hvg_genes, hvg_peaks, ratio
            )
            
            # Clean up adata_processed before saving
            # Remove large intermediate results
            for key in ['pearson_residuals_df', 'raw_pearson_residuals']:
                if key in adata_processed.uns:
                    del adata_processed.uns[key]
            
            # Save only essential information
            minimal_adata = ad.AnnData(
                X=adata_processed.X,
                obs=adata_processed.obs,
                var=adata_processed.var,
                uns={'Binary_LSI_cluster_colors': adata_processed.uns.get('Binary_LSI_cluster_colors', None)},
                obsm={'X_umap': adata_processed.obsm['X_umap']}
            )
            
            print(f"Saving processed data for ratio {ratio:.1f}")
            minimal_adata.write_h5ad(iteration_file, compression='gzip')
        
        # Save iteration-specific results
        pd.DataFrame(cell_type_stats).T.to_csv(
            os.path.join(RESULTS_DIR, f'cell_type_stats_{ratio:.1f}.csv')
        )
        pd.DataFrame(cluster_stats).T.to_csv(
            os.path.join(RESULTS_DIR, f'cluster_stats_{ratio:.1f}.csv')
        )
        
        # Create iteration-specific plots
        plot_clustering_metrics(
            pd.DataFrame(cell_type_stats).T,
            ratio,
            PLOT_DIR
        )
        
        # Store for overall analysis
        overall_results.append({
            'ratio': ratio,
            'cell_type_stats': cell_type_stats,
            'cluster_stats': cluster_stats
        })
    
    # Save overall results
    print("\nSaving overall results...")
    overall_df = pd.DataFrame(overall_results)
    overall_df.to_pickle(os.path.join(RESULTS_DIR, 'overall_results.pkl'))
    
    print(f"Analysis complete at {datetime.now():%Y-%m-%d %H:%M:%S}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        sys.exit(1)