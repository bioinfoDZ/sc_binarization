library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(ComplexHeatmap)
library(grid)


set.seed(2022)



### BINARIZATION OF DATA FROM BARYAWNO ET AL 2019
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6570562/



##########################################  define functions ############################################ 

# define function for cluster celltype mapping
celltype_cluster_map <- function(celltype_cluster, # data.frame with celltype (level) and cluster (level)
                                 cell_min_num = 10, # number of cells to be considered in the celltype
                                 singlect_thres = 0.9, # cutoff prop of single-celltype cluster
                                 multict_lothres = 0.5 # cutoff prop of multiple celltype clusters
){
  
  
  
  
  #loop thru clusters and get prop of cells
  clusts <- levels(celltype_cluster[,2])
  ctv_l <- lapply(clusts, function(cl){
    
    #get this cluster metadata
    ccc <- celltype_cluster[celltype_cluster[,2] == cl,]
    
    #get table vector
    ctv <- table(ccc[,1])
    
    #sort, remove zeros
    ctv <- ctv[ctv>0]
    ctv <- sort(ctv,decreasing = T)
    
    #cluster must have at least 10 cells?
    # celltype must contribute at lest cell_min_num cells
    ctv_cut <- ctv[ctv>cell_min_num]
    
    #get prop
    ctv <- ctv/sum(ctv)
    
    #filter pfull props by min cell cutoff
    ctv <- ctv[names(ctv) %in% names(ctv_cut)]
    
    #return ctv so we can see how it looks...
    
    return(ctv)
    
    
  })
  names(ctv_l) <- paste0('Cluster_', clusts)
  
  #based on prop of cells, try to assign celltype
  # if first is above 0.9 (singlect_thres), assign 
  # if first is below 0.9, try to assign any with more than 0.1 (multict_lothres)
  ctmapping <- sapply(clusts, function(cl){
    clustlab <- paste0('Cluster_', cl)
    
    ctv <- ctv_l[[clustlab]]
    
    #if above single cutoff thres, just pick it as celltype
    if(ctv[1] >= singlect_thres){
      cm <- names(ctv)[1]
    }
    
    #if below thres, use a combination of celltypes
    
    if(ctv[1] < singlect_thres){
      
      ctv <- ctv[ctv>multict_lothres]
      cm <- names(ctv)
      cm <- paste(cm, collapse = '_')
      
      #sometimes, no cell types pass 50%, just call it mixed
      if( cm == '' ){cm = 'Mixed'}
    }
    
    cm <- paste0(clustlab, '--', cm)
    
    
    
  })
  
  return(ctmapping)
  
}



get_top_celltype <- function(celltype_cluster,
                             cell_min_num = 10){
  
  #loop thru clusters and get prop of cells
  clusts <- levels(celltype_cluster[,2])
  cmax_v <- sapply(clusts, function(cl){
    
    #get this cluster metadata
    ccc <- celltype_cluster[celltype_cluster[,2] == cl,]
    
    #get table vector
    ctv <- table(ccc[,1])
    
    #sort, remove zeros
    ctv <- ctv[ctv>0]
    ctv <- sort(ctv,decreasing = T)
    
    #cluster must have at least 10 cells?
    # celltype must contribute at lest cell_min_num cells
    ctv_cut <- ctv[ctv>cell_min_num]
    
    #get prop
    ctv <- ctv/sum(ctv)
    
    #filter pfull props by min cell cutoff
    ctv <- ctv[names(ctv) %in% names(ctv_cut)]
    
    
    
    cmax = names(ctv[which.max(ctv)])
    
    #return ctv so we can see how it looks...
    
    return(cmax)
    
    
  })
  
  
  return(cmax_v)
}



############################################ ############################################ ############################################ 







#### first, read in full atlas, cluster, and plot celltypes... ####

atlasfile <- '~/Dropbox/data/deyou/classifier/data/bmstroma/parsed_sobj.rds'
sobjfull <- readRDS(atlasfile)



#run SCT, UMAP, clustering
sobjfull <- SCTransform(sobjfull)
sobjfull <- RunPCA(sobjfull)

ElbowPlot(sobjfull, ndims = 30)

sobjfull <- FindNeighbors(sobjfull, reduction = "pca", dims = c(1:20) )
sobjfull <- RunUMAP(sobjfull, reduction = "pca", dims = c(1:20) )
sobjfull <- FindClusters(sobjfull, resolution = 0.5)


head(sobjfull)

DimPlot(sobjfull, group.by = "seurat_clusters", label = T, repel = T) / DimPlot(sobjfull, group.by = "Celltype", label = T, repel = T)


#get cluster-celltype map
celltype_cluster <- sobjfull@meta.data[,c('Celltype', 'seurat_clusters')]
cellmap_full <- celltype_cluster_map(celltype_cluster)

cellmap_full



sobjfull$cluster_celltype <- plyr::mapvalues(sobjfull$seurat_clusters, 
                                             from = levels(sobjfull$seurat_clusters),
                                             to = cellmap_full
)


#sort cellmap by underlying celltype
cellmap_full_ctsort <- cellmap_full[str_order(str_split_fixed(cellmap_full,'--', n = 2)[,2], numeric = T)]
if(any(grepl('Mixed', cellmap_full_ctsort))){
  mixedones <- cellmap_full_ctsort[grepl('Mixed', cellmap_full_ctsort)]
  mixedones <- str_sort(mixedones, numeric = T)
  wom <- cellmap_full_ctsort[!grepl('Mixed', cellmap_full_ctsort)]
  cellmap_full_ctsort <- c(wom, mixedones)
}
sobjfull$cluster_celltype_ctsort <- factor(sobjfull$cluster_celltype, levels = cellmap_full_ctsort)



DimPlot(sobjfull, group.by = 'cluster_celltype', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()
DimPlot(sobjfull, group.by = 'cluster_celltype_ctsort', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()


beepr::beep()



#get the cmax sorted levels cluster factor
# sort clusters by top celltype...
# celltype_cluster <- sobjfull@meta.data[,c('Celltype', 'seurat_clusters')]
topcts <- get_top_celltype(celltype_cluster)
named_clusts <- levels(sobjfull$cluster_celltype)

ctdf <- data.frame(cluster = names(topcts), celltype = topcts, named_clusts = named_clusts)
ctdf <- ctdf[str_order(ctdf$celltype, numeric = T),]


sobjfull$cluster_celltype_cmaxsort <- factor(sobjfull$cluster_celltype, levels = ctdf$named_clusts)
DimPlot(sobjfull, group.by = 'cluster_celltype_cmaxsort', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()



# 
# cellmap_full_ctsort <- cellmap_full[str_order(str_split_fixed(cellmap_full,'--', n = 2)[,2], numeric = T)]
# if(any(grepl('Mixed', cellmap_full_ctsort))){
#   mixedones <- cellmap_full_ctsort[grepl('Mixed', cellmap_full_ctsort)]
#   mixedones <- str_sort(mixedones, numeric = T)
#   wom <- cellmap_full_ctsort[!grepl('Mixed', cellmap_full_ctsort)]
#   cellmap_full_ctsort <- c(wom, mixedones)
# }








### binarize matrix; 0s or 1s ###
# do this with a mask object


mat <- sobjfull@assays$RNA@counts

mask <- (mat > 0)
mat[mask] <- 1










### try the normal pipeline after binarizing ###


sobj <- CreateSeuratObject(mat, project = 'atlas')
sobj$Celltype <- sobjfull$Celltype
sobj$Full_SCT_snn_res.0.5 <- sobjfull$SCT_snn_res.0.5


#pipeline:
# var features --> use ALL genes
# pca
# graph
# cluster
# umap...



#instead of finding var features, use all genes
# update Nov 2023 - we still use var features now
sobj <- FindVariableFeatures(sobj)
# sobj@assays$RNA@var.features <- rownames(sobj) #no longer works in Seur v5
# VariableFeatures(sobj) <- rownames(sobj)

# Seurat PCA wants scaledata; coerce...
# NO NEED TO SCALE WITH JUST 1 and 0
# get the var featues and densify matrix?
# mat <- sobj@assays$RNA@counts #no longer works in Seur v5
mat <- GetAssayData(sobj, assay = 'RNA', layer = 'counts')
mat <- mat[match(VariableFeatures(sobj), rownames(mat)),] #doesn't do anything if using all features as var features


### PCA: may take a while ###



## in Seurat we are forced to input a dense matrix for PCA... ugh
# that was in Seur v4
# sobj@assays$RNA@scale.data <- as.matrix(mat)  #no longer works in Seur v5

#in Seur v5 sparse mat? --> nope still needs to be full matrix...
# sobj@assays$RNA$scale.data <- as.matrix(mat)
# 
# sobj <- Seurat::RunPCA(object = sobj, verbose = F,)
# 
# #slim it down
# sobj <- DietSeurat(sobj, counts = F, data = T, scale.data = F, dimreducs = 'pca')


# try own pca...
# pca <- prcomp(t(mat))
# also very slow...


## try using irlba

pca <- irlba::irlba( t(mat), nv=50)

#get the "PCs", ie cell embeddings, as in prcomp_irlba retx
embs <- pca$u %*% diag(pca$d)
colnames(embs) <- paste0('PC', 1:ncol(embs))
rownames(embs) <- colnames(mat)

#get standard dev
sdev <- pca$d / sqrt(ncol(mat) - 1)



# pca_prcomp_irlba <- irlba::prcomp_irlba(t(mat), n=50, center = F, scale. = F, verbose=T)
# takes so much longer than irlba, not sure why

sobj[['pca']] <- CreateDimReducObject(embeddings = embs, key = 'pca', stdev = sdev)


rm(pca,embs,sdev,mat, mask); gc(full=T)



#test against seurat PCA

# sobj2 <- sobj
# sobj2@assays$RNA$scale.data <- as.matrix(mat)
# sobj2 <- Seurat::RunPCA(object = sobj2, verbose = F,)
# beepr::beep()
# 
# (ElbowPlot(sobj) + ggtitle('irlba')) + (ElbowPlot(sobj2) + ggtitle('Seurat RunPCA'))
# DimPlot(sobj) + DimPlot(sobj2)
# perfect match, though signs are flipped, doesn't matter



#pick number of PCs
ElbowPlot(sobj, ndims = 50) + ggtitle('irlba')

dims <- c(1:15)

sobj <- Seurat::FindNeighbors(object = sobj, dims = dims, verbose = F)
sobj <- RunUMAP(sobj, dims = dims)

sobj <- Seurat::FindClusters(object = sobj, resolution = 0.5, verbose = F, algorithm = 1)


head(sobj)

# DimPlot(sobj, group.by = "seurat_clusters", label = T, repel = T) / (DimPlot(sobj, group.by = "Celltype", label = T, repel = T) + DimPlot(sobj, group.by = "Full_SCT_snn_res.0.5", label = T, repel = T))


#get cluster-celltype map
celltype_cluster <- sobj@meta.data[,c('Celltype', 'seurat_clusters')]
cellmap_bin <- celltype_cluster_map(celltype_cluster)

cellmap_full
cellmap_bin

sobj$cluster_celltype <- plyr::mapvalues(sobj$seurat_clusters, 
                                         from = levels(sobj$seurat_clusters),
                                         to = cellmap_bin
)



#sort cellmap by underlying celltype
cellmap_full_ctsort <- cellmap_bin[str_order(str_split_fixed(cellmap_bin,'--', n = 2)[,2], numeric = T)]
if(any(grepl('Mixed', cellmap_full_ctsort))){
  mixedones <- cellmap_full_ctsort[grepl('Mixed', cellmap_full_ctsort)]
  mixedones <- str_sort(mixedones, numeric = T)
  wom <- cellmap_full_ctsort[!grepl('Mixed', cellmap_full_ctsort)]
  cellmap_full_ctsort <- c(wom, mixedones)
}
sobj$cluster_celltype_ctsort <- factor(sobj$cluster_celltype, levels = cellmap_full_ctsort)



DimPlot(sobj, group.by = 'cluster_celltype', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()
DimPlot(sobj, group.by = 'cluster_celltype_ctsort', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()


beepr::beep()



#get the cmax sorted levels cluster factor
# sort clusters by top celltype...
# celltype_cluster <- sobjfull@meta.data[,c('Celltype', 'seurat_clusters')]
topcts <- get_top_celltype(celltype_cluster)
named_clusts <- levels(sobj$cluster_celltype)

ctdf <- data.frame(cluster = names(topcts), celltype = topcts, named_clusts = named_clusts)
ctdf <- ctdf[str_order(ctdf$celltype, numeric = T),]


sobj$cluster_celltype_cmaxsort <- factor(sobj$cluster_celltype, levels = ctdf$named_clusts)
DimPlot(sobj, group.by = 'cluster_celltype_cmaxsort', reduction = 'umap', label = T, repel = T) + SDAP::theme_dimplot()















### try LSI...?
# https://www.science.org/doi/10.1126/science.aab1601
# it's in the Signac package...

mat <- GetAssayData(sobj, assay = 'RNA', layer = 'counts')
tfidf <- Signac::RunTFIDF(mat)
lsi <- Signac::RunSVD(tfidf)




sobjlsi <- sobj
sobjlsi[['lsi']] <- lsi

# 
# DimPlot(sobjlsi, reduction = 'lsi', group.by = 'Celltype', label = T, repel = T)
# DimPlot(sobjlsi, reduction = 'pca', group.by = 'Celltype', label = T, repel = T)


Signac::DepthCor(sobjlsi)
Seurat::ElbowPlot(sobjlsi)

#exclude PC1?



sobjlsi <- RunUMAP(sobjlsi, reduction = 'lsi', dims=2:10, reduction.name = 'umap_lsi', reduction.key = 'umap_lsi')
sobjlsi <- RunUMAP(sobjlsi, reduction = 'lsi', dims=1:10, reduction.name = 'umap_lsi_wdim1', reduction.key = 'umap_lsi_wdim1')




DimPlot(sobjlsi, reduction = 'umap', group.by = 'Celltype', label = T, repel = T)
DimPlot(sobjlsi, reduction = 'umap_lsi', group.by = 'Celltype', label = T, repel = T)
DimPlot(sobjlsi, reduction = 'umap_lsi_wdim1', group.by = 'Celltype', label = T, repel = T)


#proceed w/o dim1 for now


sobjlsi <- FindNeighbors(object = sobjlsi, reduction = 'lsi', dims = 2:10)
# sobjlsi <- FindClusters(object = sobjlsi)
sobjlsi <- FindClusters(object = sobjlsi, resolution = 0.5)

head(sobjlsi)

#get cluster-celltype map
celltype_cluster <- sobjlsi@meta.data[,c('Celltype', 'RNA_snn_res.0.5')]
cellmap_lsi <- celltype_cluster_map(celltype_cluster)

cellmap_full
cellmap_bin
cellmap_lsi

sobjlsi$cluster_celltype <- plyr::mapvalues(sobjlsi$seurat_clusters, 
                                            from = levels(sobjlsi$seurat_clusters),
                                            to = cellmap_lsi
)


sobjlsi$cluster_res0.5_celltype <- sobjlsi$cluster_celltype

#sort cellmap by underlying celltype
cellmap_full_ctsort <- cellmap_lsi[str_order(str_split_fixed(cellmap_lsi,'--', n = 2)[,2], numeric = T)]
if(any(grepl('Mixed', cellmap_full_ctsort))){
  mixedones <- cellmap_full_ctsort[grepl('Mixed', cellmap_full_ctsort)]
  mixedones <- str_sort(mixedones, numeric = T)
  wom <- cellmap_full_ctsort[!grepl('Mixed', cellmap_full_ctsort)]
  cellmap_full_ctsort <- c(wom, mixedones)
}
sobjlsi$cluster_celltype_ctsort <- factor(sobjlsi$cluster_celltype, levels = cellmap_full_ctsort)



DimPlot(sobjlsi, group.by = 'cluster_celltype', reduction = 'umap_lsi', label = T, repel = T) + SDAP::theme_dimplot()
DimPlot(sobjlsi, group.by = 'cluster_celltype_ctsort', reduction = 'umap_lsi', label = T, repel = T) + SDAP::theme_dimplot()

beepr::beep()



#get the cmax sorted levels cluster factor
# sort clusters by top celltype...
# celltype_cluster <- sobjlsi@meta.data[,c('Celltype', 'seurat_clusters')]
topcts <- get_top_celltype(celltype_cluster)
named_clusts <- levels(sobjlsi$cluster_celltype)

ctdf <- data.frame(cluster = names(topcts), celltype = topcts, named_clusts = named_clusts)
ctdf <- ctdf[str_order(ctdf$celltype, numeric = T),]


sobjlsi$cluster_celltype_cmaxsort <- factor(sobjlsi$cluster_celltype, levels = ctdf$named_clusts)
DimPlot(sobjlsi, group.by = 'cluster_celltype_cmaxsort', reduction = 'umap_lsi', label = T, repel = T) + SDAP::theme_dimplot()













#save all outs
saveRDS(sobjfull, 'outs/boneatlas/sobjfull.rds')
saveRDS(sobj, 'outs/boneatlas/sobj_binary.rds')
saveRDS(sobjlsi, 'outs/boneatlas/sobj_tfidf_lsi.rds')


cellmaplist <- list(sobjfull_SCT_snn_res.0.5 = cellmap_full,
                    sobj_RNA_snn_res.0.5 = cellmap_bin,
                    sobjlsi_RNA_snn_res.0.5 = cellmap_lsi
)
saveRDS(cellmaplist, 'outs/cellmaplist.rds')




## READ IN AND PARSE PAPER MARKERS ##
pubm <- as.data.frame(readxl::read_excel('outs/boneatlas/NIHMS1529101-supplement-8.xlsx'))
colnames(pubm) <- pubm[1,] ; pubm <- pubm[-1,]


cluster_celltype <- sobjfull@meta.data[,c('Celltype', 'Cluster')]
cts <- get_top_celltype(cluster_celltype)

#get celltype from clusters (paper)
pubm$celltype <- pubm$cluster
pubm$celltype <- factor(pubm$celltype, levels = names(cts))
pubm$celltype <- plyr::mapvalues(pubm$celltype, from = names(cts), cts)

#numeric...
pubm$p_val <- as.numeric(pubm$p_val)
pubm$avg_logFC <- as.numeric(pubm$avg_logFC)
pubm$pct.1 <- as.numeric(pubm$pct.1)
pubm$pct.2 <- as.numeric(pubm$pct.2)
pubm$p_val_adj <- as.numeric(pubm$p_val_adj)




saveRDS(pubm, 'outs/boneatlas/PARSED_MARKERS_NIHMS1529101-supplement-8.rds')


beepr::beep()




