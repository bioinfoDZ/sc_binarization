import numpy as np
import scanpy as sc
import multigrate as mtg
import muon as mu
import os
import csv
import pandas as pd
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, mutual_info_score, fowlkes_mallows_score, silhouette_score
from anndata import AnnData

mdata = mu.read('index_1_sigma_0.8_diff_0.3_noise_0.5_processed_muon_with_mofa.h5mu')

mu.pp.intersect_obs(mdata)

rna = mdata['rna']
atac = mdata['ATAC']

rna = rna[:,rna.var.query('highly_variable == True').index]

print(rna)

atac = atac[:,atac.var.query('highly_variable == True').index]

print(atac)


adata = mtg.data.organize_multiome_anndatas(
    adatas = [[rna], [atac]],           # a list of anndata objects per modality, RNA-seq always goes first
    layers = [[None], [None]],    # if need to use data from .layers, if None use .X
)
adata.obs['pop'] = rna.obs['pop']

mtg.model.MultiVAE.setup_anndata(
    adata
)

model = mtg.model.MultiVAE(
    adata,
    losses=['mse', 'mse'],
)

model.train(max_epochs=400)

#model.plot_losses()

model.get_latent_representation()
adata

sc.pp.neighbors(adata, use_rep='latent')
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5, flavor="igraph" , n_iterations=2)

sc.pl.umap(adata, color=['pop', 'leiden'], frameon=False, ncols=1, save='sim_index1_umaps.png')

adata.write('sim_multigrate.h5ad')
adata.obs[['pop','leiden']].to_csv('multigrate_obs.csv')

model.plot_losses(save='sim_loss.png')




def calculate_metrics_and_write_to_file(data_list, output_file='metrics_a.csv'):
    def cluster_accuracy(adata, test_category, ground_truth):
        crosstab = pd.crosstab(adata.obs[test_category], adata.obs[ground_truth])
        accuracy = crosstab.max(axis=1) / crosstab.sum(axis=1)
        return accuracy.mean()

    # Data to write to the CSV file
    metrics_data = []

    for data in data_list:
        adata = data['adata']
        col1 = data['col1']
        col2 = data['col2']
        label = data['label']

        ari = adjusted_rand_score(adata.obs[col1], adata.obs[col2])
        ami = adjusted_mutual_info_score(adata.obs[col1], adata.obs[col2])
        mi = mutual_info_score(adata.obs[col1], adata.obs[col2])
        fm = fowlkes_mallows_score(adata.obs[col1], adata.obs[col2])
        acc = cluster_accuracy(adata, col2, col1)
        
        # Use precomputed PCA results for silhouette score calculation
        if 'X_pca' in adata.obsm:
            adata_pca = adata.obsm['X_pca']
        else:
            #raise ValueError(f"AnnData object {label} does not contain PCA results in 'obsm['X_pca']'")
            sc.pp.pca(adata)
            adata_pca = adata.obsm['X_pca']

        silhouette = silhouette_score(adata_pca, adata.obs[col2], metric='euclidean')

        metrics_data.append([label, ari, ami, mi, fm, acc, silhouette])

    # Write metrics to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Label', 'ARI', 'AMI', 'MI', 'FM', 'Accuracy', 'Silhouette'])  # Header
        writer.writerows(metrics_data)

data_list = [
   {'adata': adata, 'col1': 'pop', 'col2': 'leiden', 'label': 'index_1'}]

calculate_metrics_and_write_to_file(data_list)

