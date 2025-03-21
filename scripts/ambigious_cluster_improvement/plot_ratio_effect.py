# plot_ratio_effect_atac.py
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
DATA_DIR = Path("atac_outputs/processed_data")  # RNA data only
METADATA_PATH = "input_data/pbmc10k_celltypes.csv"
AMBIGUOUS_TYPES = ["intermediate mono", "CD14 mono"]
PERCENTAGES = np.arange(0, 1.1, 0.1)
CLUSTER_KEY = 'Binary_LSI_cluster'
STANDARD = 'atac:celltype'  # Single standard

def get_non_ambiguous_types(metadata):
    """Dynamically identify non-ambiguous RNA cell types"""
    rna_types = set(metadata['rna:celltype'].dropna())
    return list(rna_types - set(AMBIGUOUS_TYPES))

def calculate_purity(adata, metadata):
    """Calculate purity against RNA cell types only"""
    cluster_stats = {}

    if CLUSTER_KEY not in adata.obs:
        raise KeyError(f"Missing cluster key: {CLUSTER_KEY}")

    clusters = adata.obs[CLUSTER_KEY].astype('category').cat.categories

    for cluster in clusters:
        cluster_mask = adata.obs[CLUSTER_KEY] == cluster
        cell_ids = adata.obs[cluster_mask].index.intersection(metadata.index)

        true_labels = metadata.loc[cell_ids, STANDARD].dropna()
        if true_labels.empty:
            continue

        dominant_type = true_labels.mode()[0]
        purity = (true_labels == dominant_type).mean()

        cluster_stats[cluster] = {
            'dominant_type': dominant_type,
            'purity': purity,
            'size': len(true_labels)
        }

    purity_data = {}
    non_ambiguous = get_non_ambiguous_types(metadata)

    for ct in non_ambiguous:
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

def plot_results(results, output_path):
    """Simplified visualization for RNA-only analysis"""
    plot_df = pd.DataFrame(results)
    non_ambiguous = get_non_ambiguous_types(pd.read_csv(METADATA_PATH))

    plt.figure(figsize=(12, 6))
    ax = plt.gca()
    colors = plt.cm.tab20.colors

    for i, ct in enumerate(non_ambiguous):
        subset = plot_df[plot_df['Cell Type'] == ct]
        ax.plot(
            subset['Percentage'],
            subset['Purity'],
            color=colors[i % 20],
            marker='o',
            linestyle='-',
            linewidth=2,
            markersize=8,
            label=ct
        )

    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=False
    )

    ax.set_title("ATAC+RNA(ATAC Standard)", pad=15)
    ax.set_xlabel("Percentage of Genes Added", labelpad=12)
    ax.set_ylabel("Congregated Score", labelpad=12)
    ax.set_xticks(np.arange(0, 110, 10))
    #ax.set_ylim(0.7, 1.05)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

def main():
    metadata = pd.read_csv(METADATA_PATH, index_col=0)
    results = []

    for pct in PERCENTAGES:
        pct_str = f"{pct:.1f}"
        try:
            adata_path = DATA_DIR / f"iteration_{pct_str}.h5ad"
            adata = sc.read_h5ad(adata_path)

            purity = calculate_purity(adata, metadata)

            for ct, value in purity.items():
                results.append({
                    'Percentage': pct*100,
                    'Cell Type': ct,
                    'Purity': value
                })

        except Exception as e:
            print(f"Error processing {pct_str}: {str(e)}")
            continue

    plot_results(results, "subfigure6.png")

if __name__ == "__main__":
    main()
