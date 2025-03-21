import muon as mu
import mudata
import scanpy as sc
import numpy as np
import torch
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import scmomat
import scmomat.model as model
import scmomat.utils as utils
import scmomat.umap_batch as umap_batch
from umap import UMAP
import scipy.sparse as sp
import pandas as pd
import os
from sklearn import metrics
import csv
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, mutual_info_score, fowlkes_mallows_score, silhouette_score
from anndata import AnnData
import numpy as np




# Set device for model computations
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Path to the muon file
muon_file_path = 'index_1_sigma_0.8_diff_0.3_noise_0.5_processed_muon_with_mofa.h5mu'

muon_file = mu.read(muon_file_path)
mu.pp.intersect_obs(muon_file)
rna_muon = muon_file['rna']
atac_muon = muon_file['ATAC']
print(rna_muon.var['highly_variable'].sum())
print(len(rna_muon.obs.index))

rna_hvg = rna_muon.var[rna_muon.var['highly_variable']].index
atac_hv = atac_muon.var[atac_muon.var['highly_variable']].index



data_dir = "index_1_sigma_0.8_diff_0.3_noise_0.5"
result_dir = "results_sim_index/"
os.makedirs(result_dir, exist_ok=True)


adata = mu.read_10x_mtx(data_dir)
adata.var_names_make_unique()
adata
mu.pp.intersect_obs(adata)
adata


rna_adata = adata['rna']

rna_counts = rna_adata.X  
# Properly convert sparse matrix to dense array if necessary
rna_counts = rna_counts.toarray() if sp.issparse(rna_counts) else rna_counts
rna_processed = scmomat.preprocess(rna_counts, modality="RNA", log=False)

atac_data = adata['ATAC']

atac_counts = atac_data.X
# Convert sparse matrix to dense array if necessary
atac_counts = atac_counts.toarray() if sp.issparse(atac_counts) else atac_counts

atac_processed = scmomat.preprocess(atac_counts, modality="ATAC")

feats_name = {"rna": rna_adata.var.index, "atac": atac_data.var.index}


# Wrap the processed RNA and ATAC data into a list
rna_processed_batches = [rna_processed]  # Even if it's just one batch, it needs to be a list
atac_processed_batches = [atac_processed]


try:
    scmomat_model = model.scmomat_model(counts={
        'feats_name': feats_name,
        'nbatches':1,
        'rna': rna_processed_batches,
        'atac': atac_processed_batches
    }, K=30, device=device)
   
    T = 4000
    losses = scmomat_model.train_func(T=T)


    # Plot the training losses
    plt.figure(figsize=(10, 5))
    plt.plot(losses, label='Training Loss')
    plt.title('Training Loss Over Time')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(os.path.join(result_dir, 'sim_training_loss.png'))
    plt.show()
    
    # Define parameters for post-processing
    n_neighbors = 30
    njobs = 8
    

    # Execute post-processing to obtain knn indices and distances
    
    alt_zs = scmomat_model.extract_cell_factors()

    print(type(alt_zs))
    print(len(alt_zs))
    print(type(alt_zs[0]))
    print(alt_zs[0].shape)

    print('start post graph calculation')
    knn_indices, knn_dists = scmomat.calc_post_graph(alt_zs, n_neighbors, njobs=njobs)
    print('post graph calculated')
    umap_op = UMAP(n_components = 2, n_neighbors = 15, min_dist = 0.2, random_state = 0)

    x_umap = umap_op.fit_transform(np.concatenate(alt_zs, axis=0))
    print(x_umap.shape)
    rna_muon.obsm['X_momat_umap'] = x_umap

    embeddings = np.concatenate(alt_zs, axis=0)
    pca = PCA(n_components=15)  # Adjust the number of components as needed
    embeddings_pca = pca.fit_transform(embeddings)


    print("umap calculated")

    # Clustering using the Leiden algorithm
    resolution = 0.6
    labels_leiden = scmomat.leiden_cluster(X=None, knn_indices=knn_indices, knn_dists=knn_dists, resolution=resolution)
    print(labels_leiden)
    print(type(labels_leiden))
    

    x_umap_post = scmomat.calc_umap_embedding(knn_indices = knn_indices, knn_dists = knn_dists, n_components = 2, n_neighbors = 15, min_dist = 0.20, random_state = 0)
    rna_muon.obsm['X_momat_post_umap'] = x_umap_post
    print(x_umap_post.shape)
    

    def calculate_metrics_and_write_to_file(data_list, output_file='metrics_momat_precomputed.csv'):
        metrics_data = []

        def cluster_accuracy(test_category, ground_truth):
            crosstab = pd.crosstab(test_category, ground_truth)
            accuracy = crosstab.max(axis=1) / crosstab.sum(axis=1)
            return accuracy.mean()


        for data in data_list:
            col1 = data['col1']
            col2 = data['col2']
            label = data['label']

        ari = adjusted_rand_score(col1, col2)
        ami = adjusted_mutual_info_score(col1, col2)
        mi = mutual_info_score(col1, col2)
        fm = fowlkes_mallows_score(col1, col2)
        acc = cluster_accuracy(col1, col2)
        silhouette = silhouette_score(embeddings_pca, labels_leiden, metric='euclidean')

        metrics_data.append([label, ari, ami, mi, fm, silhouette, acc])

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Label', 'ARI', 'AMI', 'MI', 'FM', 'Silhouette','Accuracy'])  # Header
            writer.writerows(metrics_data)

    data_list = [{'col1': rna_muon.obs['pop'], 'col2': labels_leiden, 'label': 'index_1'}]

    calculate_metrics_and_write_to_file(data_list)

    # Plot the UMAP results
    scmomat.plot_latent(x_umap, annos=labels_leiden, mode="joint", axis_label="UMAP", markerscale=6, s=5, label_inplace=True, alpha=0.7,save='sim_momat_umap_leiden.png')
    scmomat.plot_latent(x_umap, annos=rna_muon.obs["pop"], mode="joint", axis_label="UMAP", markerscale=6, s=5, label_inplace=True, alpha=0.7,save='sim_momat_umap_rna.png')

    scmomat.plot_latent(x_umap_post, annos = labels_leiden, mode = "joint", save = 'index_sim_post_umap_leiden.png',\
                    axis_label = "UMAP", markerscale = 6, s = 5, label_inplace = False, alpha = 0.7)

    scmomat.plot_latent(x_umap_post, annos = rna_muon.obs["pop"], mode = "joint", save = 'index_sim_post_umap_celltype.png',\
                    axis_label = "UMAP", markerscale = 6, s = 5, label_inplace = False, alpha = 0.7)
    

    rna_muon.obs['labels_leiden'] = list(labels_leiden)
    rna_muon.obs.to_csv('momat_obs.cv')

    rna_muon.obs['labels_leiden'] = rna_muon.obs['labels_leiden']

    sc.pl.embedding(rna_muon, color='pop', basis='X_momat_umap', save='index_sim_scanpy_umap.png')
    sc.pl.embedding(rna_muon, color='pop', basis='X_momat_post_umap', save='index_sim_scanpy_post_umap.png')


    sc.pl.embedding(rna_muon, color=labels_leiden, basis='X_momat_umap', save='index_sim_scanpy_umap_leiden.png')
    sc.pl.embedding(rna_muon, color=labels_leiden, basis='X_momat_post_umap', save='index_sim_scanpy_post_umap_leiden.png')

except Exception as e:
    print("Error in model processing:", e)

print("scMoMaT analysis completed or failed with error.")

