import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
#warnings.filterwarnings('ignore')

# Configuration
DATA_DIR = Path("atac_outputs/processed_data")
METADATA_PATH = "input_data/pbmc10k_celltypes.csv"
AMBIGUOUS_TYPES = ["intermediate mono", "CD14 mono"]
PERCENTAGES = np.arange(0, 1.1, 0.1)
STANDARDS = ['rna:celltype', 'atac:celltype']
CLUSTER_KEY = 'Binary_LSI_cluster'  # Verified from analyze_peaks.py

def calculate_dual_purity(adata, metadata, standard):
    """Calculate purity against specified annotation standard with error handling"""
    try:
        cluster_stats = {}

        # Verify cluster key exists
        if CLUSTER_KEY not in adata.obs:
            raise KeyError(f"Cluster key '{CLUSTER_KEY}' not found in adata.obs")

        # Get unique clusters
        clusters = adata.obs[CLUSTER_KEY].astype('category').cat.categories

        for cluster in clusters:
            cluster_mask = adata.obs[CLUSTER_KEY] == cluster
            cluster_cells = adata.obs[cluster_mask]
            cell_ids = cluster_cells.index.intersection(metadata.index)

            if len(cell_ids) == 0:
                print(f"Skipping cluster {cluster} - no cells match metadata")
                continue

            # Get ground truth from metadata
            true_labels = metadata.loc[cell_ids, standard].dropna()

            if len(true_labels) == 0:
                print(f"Skipping cluster {cluster} - no labels for {standard}")
                continue

            dominant_type = true_labels.value_counts().idxmax()
            purity = (true_labels == dominant_type).mean()

            cluster_stats[cluster] = {
                'dominant_type': dominant_type,
                'purity': purity,
                'size': len(cell_ids)
            }

        # Calculate weighted averages for ambiguous types
        purity_data = {}
        for ct in AMBIGUOUS_TYPES:
            relevant_clusters = [
                c for c, stats in cluster_stats.items()
                if stats['dominant_type'] == ct
            ]

            if relevant_clusters:
                total_size = sum(cluster_stats[c]['size'] for c in relevant_clusters)
                weighted_purity = sum(
                    cluster_stats[c]['purity'] * cluster_stats[c]['size']
                    for c in relevant_clusters
                ) / total_size
            else:
                weighted_purity = 0.0

            purity_data[ct] = weighted_purity

        return purity_data

    except Exception as e:
        print(f"Error in calculate_dual_purity: {str(e)}")
        return {ct: 0.0 for ct in AMBIGUOUS_TYPES}

def main():
    # Load metadata
    metadata = pd.read_csv(METADATA_PATH, index_col=0)

    # Process all percentage points
    results = []
    for pct in PERCENTAGES:
        pct_str = f"{pct:.1f}"
        try:
            adata = sc.read_h5ad(DATA_DIR / f"iteration_{pct_str}.h5ad")

            # Verify required data exists
            if 'X_umap' not in adata.obsm:
                print(f"Warning: Missing UMAP coordinates in {pct_str}")
                adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))

            # Calculate purity against both standards
            for standard in STANDARDS:
                purity = calculate_dual_purity(adata, metadata, standard)
                if purity is None:
                    raise ValueError("Purity calculation returned None")

                for ct, value in purity.items():
                    results.append({
                        'Percentage': pct*100,
                        'Cell Type': ct,
                        'Standard': standard.replace(':celltype', '').upper(),
                        'Purity': value
                    })

        except Exception as e:
            print(f"Error processing {pct_str}: {str(e)}")
            continue

    # Create plot
    plot_df = pd.DataFrame(results)

    # Visualization
    # Visualization
    plt.figure(figsize=(14, 7))
    ax = plt.gca()

    # New color and style mapping
    COLOR_MAP = {
        'intermediate mono': '#1f77b4',  # Blue
        'CD14 mono': '#ff7f0e'           # Orange
    }

    LINE_STYLES = {
        'RNA': '-',
        'ATAC': '--'
    }

    MARKERS = {
        'RNA': 'o',
        'ATAC': 's'
    }

    # Plot both standards for each cell type
    for ct in AMBIGUOUS_TYPES:
        for standard in ['RNA', 'ATAC']:
            subset = plot_df[(plot_df['Cell Type'] == ct) & 
                           (plot_df['Standard'] == standard)]

            ax.plot(
                subset['Percentage'],
                subset['Purity'],
                color=COLOR_MAP[ct],
                linestyle=LINE_STYLES[standard],
                marker=MARKERS[standard],
                markersize=8,
                markerfacecolor='white',
                markeredgewidth=1.5,
                linewidth=2,
                label=f"{ct} ({standard})"
            )

    # Create unified legend
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))  # Remove duplicates
    ax.legend(
        unique_labels.values(),
        unique_labels.keys(),
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=False,
        fontsize=10
    )

    # Add annotation for line styles
    ax.text(
        1.02, 0.4, 
        "Line Style Key:\n"
        "─ RNA Standard\n" 
        "▬▬▬ ATAC Standard",
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment='top'
    )

    ax.set_title("ATAC+RNA", pad=20)
    ax.set_xlabel("Percentage of Genes Added", labelpad=12)
    ax.set_ylabel("Congregated Score", labelpad=12)
    plt.xticks(np.arange(0, 110, 10))
    plt.ylim(0.4, 1.0)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("dual_standard_purity_comparison_rna_atac2b.pdf", 
               dpi=300, 
               bbox_inches='tight',
               transparent=True)


if __name__ == "__main__":
    main()