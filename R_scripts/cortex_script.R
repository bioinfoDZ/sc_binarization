library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(ComplexHeatmap)
library(grid)
library(Signac)


set.seed(2022)



### BINARIZATION OF DATA FROM ZHU ET AL SCI ADV 2023 
# https://pubmed.ncbi.nlm.nih.gov/37824614/








## read data
# count matrices for RNA and ATAC are from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204684

## read RNA
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204683/suppl/GSE204683%5Fcount%5Fmatrix.RDS.gz


rna <- readRDS("data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/expression/GEO_GSE204684/rna/GSE204683_count_matrix.RDS")



## read ATAC
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204682
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204682/suppl/GSE204682%5Fcount%5Fmatrix.RDS.gz



atac <- readRDS("data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/expression/GEO_GSE204684/atac/GSE204682_count_matrix.RDS")



## read in metadata from broad portal
# https://singlecell.broadinstitute.org/single_cell/study/SCP1859/multi-omic-profiling-of-the-developing-human-cerebral-cortex-at-the-single-cell-level#study-summary

md <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/metadata/scp_metadata.csv')
md <- md[-1,]

# # also read embeddings --> DO THIS LATER
# # RNA
# emb_rna <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusterrnainfo.csv')
# emb_rna <- emb_rna[-1,]
# 
# # ATAC  
# emb_atac <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusteratacinfo.csv')
# emb_atac <- emb_atac[-1,]
# 
# # WNN  
# emb_wnn <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusterwnninfo.csv')
# emb_wnn <- emb_wnn[-1,]


# embs <- data.frame(barcodes = emb_rna$NAME,
#                    RNA_UMAP_1 = emb_rna$X, RNA_UMAP_2 = emb_rna$Y,
#                    ATAC_UMAP_1 = emb_atac$X, ATAC_UMAP_2 = emb_atac$Y,
#                    WNN_UMAP_1 = emb_wnn$X, WNN_UMAP_2 = emb_wnn$Y
#                    )
# 
# rm(emb_atac,emb_rna, emb_wnn)
# gc(full = T)







### try to process as Seurat ###
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac

## PAPER ALREADY DID THIS, JUST BINARIZE AND TRY TO MATCH CLUSTER NUM


# 
# ### RNA preproc
# # rna <- as(rna[["RNA"]], Class = "Assay5")
# rna <- CreateSeuratObject(rna, assay = 'RNA')
# 
# rna <- NormalizeData(rna)
# rna <- FindVariableFeatures(rna)
# rna <- ScaleData(rna)
# rna <- RunPCA(rna)
# 
# # ElbowPlot(rna, ndims = 50)
# 
# rna <- RunUMAP(rna, dims = 1:30)
# 
# 
# 
# 
# ## ATAC preproc
# 
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "hg38"
# Annotation(pbmc.atac) <- annotations
# 
# 
# 













### BINARIZE ####

# binarize, conatenate rna+atac, treat as atac (tf/idf, etc)


### binarize rna matrix; 0s or 1s ###
# do this with a mask object
mask <- (rna > 0)
rna[mask] <- 1

rm(mask) ; gc(full=T)






### conatenate rna and atac ###
#keep track of feature names / types
ftt <- data.frame(feature = c(rownames(atac), rownames(rna)),
                  type = c( rep('atac', nrow(atac)), rep('rna', nrow(rna)))
)
mat <- rbind(atac, rna)
rm(atac, rna) ; gc(full = T)



### make seurat object ###
# try to coerce...
# chrom_assay <- CreateChromatinAssay(
#   counts = mat,
#   sep = c(":", "-")
#   # fragments = '../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
#   # min.cells = 10,
#   # min.features = 200
# )
sobj <- CreateSeuratObject(counts = mat, assay = 'peaks')


#add in metadata
sobj$sex <- md$sex
sobj$age <- md$age
sobj$age_group <- md$age_group
sobj$donor_id <- md$donor_id
sobj$celltype <- md$celltype



#add embeddings
# also read embeddings
# RNA
emb_rna <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusterrnainfo.csv')
emb_rna <- emb_rna[-1,]

rownames(emb_rna) <- emb_rna[,1]
emb_rna <- emb_rna[,-1]
colnames(emb_rna) <- paste0('Published_UMAP_RNA_', 1:2)
emb_rna[,1] <- as.numeric(emb_rna[,1])
emb_rna[,2] <- as.numeric(emb_rna[,2])
emb_rna <- as.matrix(emb_rna)
emb_rna <- CreateDimReducObject(emb_rna, key = 'Published_UMAP_RNA')

sobj[['Published_UMAP_RNA']] <- emb_rna


# ATAC  
emb_atac <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusteratacinfo.csv')
emb_atac <- emb_atac[-1,]

rownames(emb_atac) <- emb_atac[,1]
emb_atac <- emb_atac[,-1]
colnames(emb_atac) <- paste0('Published_UMAP_atac_', 1:2)
emb_atac[,1] <- as.numeric(emb_atac[,1])
emb_atac[,2] <- as.numeric(emb_atac[,2])
emb_atac <- as.matrix(emb_atac)
emb_atac <- CreateDimReducObject(emb_atac, key = 'Published_UMAP_ATAC')

sobj[['Published_UMAP_ATAC']] <- emb_atac

# WNN  
emb_wnn <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/cluster/clusterwnninfo.csv')
emb_wnn <- emb_wnn[-1,]

rownames(emb_wnn) <- emb_wnn[,1]
emb_wnn <- emb_wnn[,-1]
colnames(emb_wnn) <- paste0('Published_UMAP_wnn_', 1:2)
emb_wnn[,1] <- as.numeric(emb_wnn[,1])
emb_wnn[,2] <- as.numeric(emb_wnn[,2])
emb_wnn <- as.matrix(emb_wnn)
emb_wnn <- CreateDimReducObject(emb_wnn, key = 'Published_UMAP_WNN')

sobj[['Published_UMAP_WNN']] <- emb_wnn


rm(emb_atac, emb_rna, emb_wnn, md, mat)
gc(full = T)








#run tfidf and lsi
mat <- sobj@assays$peaks@layers$counts

# featuere selection; or use all? 
# sobj <- FindTopFeatures(sobj, min.cutoff = 'q0')
# mat <- mat[rownames(mat) %in% VariableFeatures(sobj),]

# tf idf
tfidf <- Signac::RunTFIDF(mat)
rownames(tfidf) <- rownames(mat)
colnames(tfidf) <- colnames(mat)


dimnames(tfidf) <- list(rownames(sobj), colnames(sobj) )

# lsi
lsi <- RunSVD(tfidf)
# not sure why it is not working...
# Running SVD
# Scaling cell embeddings
# Error in validObject(.Object) : 
#   invalid class “DimReduc” object: rownames must be present in 'cell.embeddings'
# In addition: Warning message:
#   No assay specified, setting assay as RNA by default. 

# try to just run it myself
# https://github.com/stuart-lab/signac/blob/8ecdde291676102bb3b503f48926c993354b5471/R/dimension_reduction.R#L70
# library(irlba)
# n = 50
# irlba.work = n * 3
# tol = 1e-05
# reduction.key = "LSI_"
# object <- tfidf
# 
# components <- irlba(A = t(x = object), nv = n, work = irlba.work, tol = tol) #takes a couple mins
# 
# feature.loadings <- components$v
# sdev <- components$d / sqrt(x = max(1, nrow(x = object) - 1))
# cell.embeddings <- components$u
# embed.mean <- apply(X = cell.embeddings, MARGIN = 2, FUN = mean)
# embed.sd <- apply(X = cell.embeddings, MARGIN = 2, FUN = sd)
# norm.embeddings <- t((t(cell.embeddings) - embed.mean) / embed.sd)
# rownames(x = feature.loadings) <- rownames(x = object)
# colnames(x = feature.loadings) <- paste0(
#   reduction.key, seq_len(length.out = n)
# )
# rownames(x = norm.embeddings) <- colnames(x = object)
# colnames(x = norm.embeddings) <- paste0(
#   reduction.key, seq_len(length.out = n)
# )
# reduction.data <- CreateDimReducObject(
#   embeddings = norm.embeddings,
#   loadings = feature.loadings,
#   # assay = 'peaks',
#   stdev = sdev,
#   key = reduction.key,
#   misc = components
# )

rm(components, feature.loadings, mat, norm.embeddings, object, tfidf, embed.mean, embed.sd, irlba.work, n , reduction.key, sdev, tol, cell.embeddings)
gc(full = T)

#add to sobj
sobj[['lsi']] <- lsi
rm(lsi); gc(full = T)


#check first "PC" to see if driven by "depth"
DepthCor(sobj)
ElbowPlot(sobj, ndims = 50, reduction = 'lsi')




### downstream analysis: clustering and umap
#umap
sobj <- RunUMAP(sobj, reduction = 'lsi', dims=2:20, reduction.name = 'umap_lsi', reduction.key = 'umap_lsi')

#graph
sobj <- FindNeighbors(object = sobj, reduction = 'lsi', dims = 2:20)

#clustering; try to get ~15 clusters
table(sobj$celltype)
length(table(sobj$celltype))
sobj <- FindClusters(object = sobj, resolution = 0.2, graph.name = 'RNA_snn')



reduction = 'umap_lsi'
DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'seurat_clusters') / 
  DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'celltype')










## save outs
saveRDS(sobj, 'outs/cerebralcortex/sobj_binary_tfidf_lsi.rds')
saveRDS(ftt, 'outs/cerebralcortex/feature_types.rds')

beepr::beep()
























##### UPDATE CHECK SOME MARKERS TO TRY TO UNDERSTAND THE MIXED CLUSTER #####


# read in binarized, get cluster ids
sobjbin <- readRDS("outs/cerebralcortex/sobj_binary_tfidf_lsi.rds")

binclusts <- sobjbin$RNA_snn_res.0.2

rm(sobjbin); gc(full = T)

## read RNA
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204683/suppl/GSE204683%5Fcount%5Fmatrix.RDS.gz


rna <- readRDS("data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/expression/GEO_GSE204684/rna/GSE204683_count_matrix.RDS")


## read in metadata from broad portal
# https://singlecell.broadinstitute.org/single_cell/study/SCP1859/multi-omic-profiling-of-the-developing-human-cerebral-cortex-at-the-single-cell-level#study-summary

md <- read.csv('data/BROAD_HUMAN_CEREBRAL_CORTEX/SCP1859/metadata/scp_metadata.csv')
md <- md[-1,]





#make seurat object and proc
sobj <- CreateSeuratObject(rna, meta.data = md)

rm(rna, md); gc(full = T)

sobj <- SCTransform(sobj)


audio::wait(10)
gc(full = T)
beepr::beep()


#use makrers from fig 1F
markers <- c('VIM', 'PAX6', 'HES5', 'EOMES', 'SATB2', 'SLC17A7', 'NEUROD2', 'GAD2', 'LHX6', 'VIP', 'OLIG1', 'SOX10', 'AQP4', 'OPALIN', 'PTPRC', 'CLDN5', 'PDGFRB', 'COL1A2')
markers[!(markers %in% rownames(sobj))]




sobj$binary_RNA_snn_res.0.2 <- binclusts


sobj$celltype_dotplot <- factor(sobj$celltype, levels = c('RG','IPC', 'EN-fetal-early', 'EN-fetal-late', 'EN', 'IN-fetal', 'IN-MGE', 'IN-CGE', 'OPC', 'Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Pericytes', 'VSMC'))
DotPlot(sobj, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 


#get top ct
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

celltype_cluster <- sobj@meta.data[,c('celltype_dotplot', 'binary_RNA_snn_res.0.2')]
cmax = get_top_celltype(celltype_cluster)

#try to match levs
levs <- levels(sobj$celltype_dotplot)
levs <- levs[levs%in%cmax]
cmax = factor(cmax, levels = levs)
cmax = cmax[order(cmax)]

sobj$binary_RNA_snn_res.0.2_cmax <- factor(sobj$binary_RNA_snn_res.0.2, levels = names(cmax))

DotPlot(sobj, features = markers, group.by = 'binary_RNA_snn_res.0.2_cmax')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 



#vln plots of qc
sobj$nCount_RNA <- as.numeric(sobj$nCount_RNA)
sobj$nFeature_RNA <- as.numeric(sobj$nFeature_RNA)
sobj$percent_mt <- as.numeric(sobj$percent_mt)

VlnPlot(sobj, c('nCount_RNA', 'nFeature_RNA', 'percent_mt'), ncol = 1)


## find markers for c10
# sobj <- SetIdent(sobj, value = 'binary_RNA_snn_res.0.2')
# m <- FindAllMarkers(sobj, only.pos = T); beepr::beep()
# m <- m[m$p_val_adj < 0.05,]
# m10 = m[m$cluster=='10',]
# write.csv(m, 'outs/cerebralcortex/markers_binaryclusters_RNA_snn_res.0.2_calc-from-nonbinary.csv', quote = F, row.names = F)

m <- read.csv('outs/cerebralcortex/markers_binaryclusters_RNA_snn_res.0.2_calc-from-nonbinary.csv')
m <- m[m$p_val_adj < 0.05,]
m10 = m[m$cluster=='10',]





## make two more dot plots
# 1. with all cells except cluster 10, using the paper markers, and using the paper cell types
md <- sobj@meta.data
md <- md[md$binary_RNA_snn_res.0.2 != 10,]
sobj_noc10 <- sobj[,rownames(md)]


DotPlot(sobj_noc10, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0.2,0,0,0),"cm")
  ) +
  ggtitle('No Cluster 10 cells')


# 2. with only cluster 10 cells, using paper markers, using paper cell types
md <- sobj@meta.data
md <- md[md$binary_RNA_snn_res.0.2 == 10,]
sobj_onlyc10 <- sobj[,rownames(md)]


DotPlot(sobj_onlyc10, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0.2,0,0,0),"cm")
  ) +
  ggtitle('Only Cluster 10 cells')



saveRDS(sobj, 'outs/cerebralcortex/sobj_nonbinary_rnaonly.rds')











### UPDATE: we will read in the object, try a range of ALL ATAC + INCREASING NUMBERS OF HVGs, and calcualte RAND INDEX ( + ARI)


sobj <- readRDS("outs/cerebralcortex/sobj_binary_tfidf_lsi.rds")
ftt <- readRDS('outs/cerebralcortex/feature_types.rds')


## HVGs...
mat <- GetAssayData(sobj, assay = 'peaks', layer = 'counts')

#get RNA features
rna_ftt <- ftt[ftt$type == 'rna',]
mat <- mat[rownames(mat) %in% rna_ftt$feature,]


# find HVGs
vf <- FindVariableFeatures(mat)
vf <- vf[order(vf$vst.variance, decreasing = T),]



# save original, which used all features
sobj_orig <- sobj




#add in amount of 
numhvgs <- c(0, 100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 10000, 15000, 30201)

timestart = proc.time()

gc(full = T)
hvg_loop_res <- lapply(numhvgs, function(nh){
  
  
  
  hvf_filename <- paste0('outs/cerebralcortex/hvf_loop/', 'rna_', nh, '.rds')
  
  
  if(file.exists(hvf_filename)){
    
    gc(full = T)
    
    message(nh, ' already done')
    
    gc(full = T)
    
    return(outlist)
    
  } else{
    
    
    
    message('##################\n',
            '##################\n',
            '##################\n',
            '\n\n',
            'Now Running: ', nh, '\n',
            '\n\n',
            '##################\n',
            '##################\n',
            '##################\n'
            
    )
    
    
    
    #get sobj
    sobj <- sobj_orig
    
    gc(full = T)
    
    #get hvgs for this iteration
    vf_sub <- vf[1:nh,]
    if(nh == 0){
      vf_sub <- vf[0,,drop=F]
    }
    
    #add in HVGs from binary RNA, make sure to add atac too
    atac_features <- ftt[ftt$type == 'atac',"feature"]
    varfeats <- c(rownames(vf_sub), atac_features)
    
    VariableFeatures(sobj) <- varfeats
    
    
    
    #run tfidf and lsi
    mat <- GetAssayData(sobj, assay = 'peaks', layer = 'counts')
    mat <- mat[rownames(mat) %in% varfeats,]
    
    # featuere selection; or use all? 
    # sobj <- FindTopFeatures(sobj, min.cutoff = 'q0')
    # mat <- mat[rownames(mat) %in% VariableFeatures(sobj),]
    
    # tf idf
    tfidf <- Signac::RunTFIDF(mat)
    
    
    # lsi
    lsi <- RunSVD(tfidf)
    
    
    sobj[['lsi']] <- lsi
    
    
    rm(tfidf, lsi, mat)
    gc(full = T)
    
    
    ### downstream analysis: clustering and umap
    #umap
    sobj <- RunUMAP(sobj, reduction = 'lsi', dims=2:20, reduction.name = 'umap_lsi', reduction.key = 'umap_lsi')
    
    #graph
    sobj <- FindNeighbors(object = sobj, reduction = 'lsi', dims = 2:20)
    
    #clustering; try to get ~15 clusters
    # table(sobj$celltype)
    # length(table(sobj$celltype))
    sobj <- FindClusters(object = sobj, resolution = 0.2, graph.name = 'RNA_snn')
    
    
    
    reduction = 'umap_lsi'
    # dp <- (DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'seurat_clusters') + SDAP::theme_dimplot()) / 
    #   (DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'celltype') + SDAP::theme_dimplot())
    # 
    
    
    #get num atac and rna features
    ftt_in <- ftt[ftt$feature %in% varfeats,]
    ftt_in$type <- factor(ftt_in$type, levels = c('atac', 'rna'))
    ftt_tab <- table(ftt_in$type)
    
    # #prep title of patchwork
    # dptitle <- paste0( 'ATAC features: ', ftt_tab[1], '; RNA features: ', ftt_tab[2])
    # 
    #try to get rand index, ari
    message('\n\nCALCULATE RAND INDEX\n\n\n')
    celltype_cluster <- sobj@meta.data[,c('celltype', 'seurat_clusters')]
    celltype_cluster$celltype <- factor(celltype_cluster$celltype)
    
    # rand <- fossil::rand.index(celltype_cluster[,1], celltype_cluster[,2])
    # adj.rand <- fossil::adj.rand.index(celltype_cluster[,1], celltype_cluster[,2])
    
    adj.rand <- round( pdfCluster::adj.rand.index(celltype_cluster[,1], celltype_cluster[,2]) , 2)
    
    
    
    
    ## also calculate "accuracy"
    twt <- as.matrix(table(celltype_cluster[,1], celltype_cluster[,2]))
    outs <- lapply(1:ncol(twt), function(i){
      x = twt[,i]
      
      (max(x) / sum(x)) * 100
    })
    outvec <- unlist(outs)
    names(outvec) <- colnames(twt)
    
    
    meanacc <- round(mean(outvec) , 2)
    minacc <- round(min(outvec) , 2)
    maxacc <- round( max(outvec) , 2)
    
    # dpsubtitle <- paste0('Adjusted Rand Index = ', adj.rand, '\n',
    #                      'Column Accuracy % [range] = ',  meanacc, ' [', minacc, ' - ', maxacc , ']')
    # 
    # 
    # 
    # dp <- dp + plot_annotation(title =  dptitle,
    #                            subtitle = dpsubtitle)
    
    
    #also get num clusts
    numclusts = length(levels(celltype_cluster[,2]))
    
    rm(mat, tfidf, lsi, vf_sub, ftt_in, varfeats, atac_features, celltype_cluster, rna_ftt, outs, twt)
    gc(full = T)
    
    
    df = data.frame(atac_features = ftt_tab[1],
                    rna_features = ftt_tab[2],
                    ari = adj.rand,
                    meanacc = meanacc,
                    minacc = minacc,
                    maxacc = maxacc,
                    numclusts = numclusts)
    
    
    message('\nSaving...')
    
    saveRDS(list(sobj = sobj,
                 outvec = outvec,
                 df = df),
            file = hvf_filename
    )
    
    
    rm(sobj)
    gc(full = T)
    
    
    return(list(outvec = outvec,
                df = df)
    )
    
  } #file loop
  
})


gc(full = T)

beepr::beep()

hourspassed <- (proc.time() - timestart)[3]/60/60
names(hourspassed) <- 'Hours'
hourspassed

names(hvg_loop_res) <- paste0('rna_', numhvgs)


saveRDS(hvg_loop_res, "outs/cerebralcortex/hvf_loop/rna_sumstats.rds")
hvg_loop_res <- readRDS("outs/cerebralcortex/hvf_loop/rna_sumstats.rds")



RNALOOP = hvg_loop_res


dfs <- lapply(1:length(hvg_loop_res), function(i){
  hvg_loop_res[[i]]$df
})
dfs <- dplyr::bind_rows(dfs)
rownames(dfs) <- NULL


outvecs = lapply(1:length(hvg_loop_res), function(i){
  thisres = hvg_loop_res[[i]]
  
  df = thisres$df
  ov = thisres$outvec
  
  clusteraccdf = data.frame(rna_features = df$rna_features,
                            cluster = names(ov),
                            clusteraccuracy = ov)
  
  clusteraccdf
  
  
})
outvecs <- dplyr::bind_rows(outvecs)
rownames(outvecs) <- NULL


# make plot of rna numbers versus mean acc, rna numbers vs ari


write.csv(dfs, 'outs/cerebralcortex/hvf_loop/RNA_DF.csv', row.names = F, quote = F)

dfm <- dfs
dfm$rna_features <- prettyNum(dfm$rna_features, big.mark = ',')
dfm$rna_features <- factor(dfm$rna_features, levels = dfm$rna_features)

rna_ari <- ggplot(dfm, aes(rna_features, ari))+
  geom_point()+
  theme_light()+
  ylab('Adjusted Rand Index') + xlab('Number of included RNA features sorted by variance')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying RNA features, effect on A.R.I.')




rna_acc <- ggplot(dfm, aes(rna_features, meanacc))+
  geom_point()+
  theme_light()+
  ylab('Mean Cluster Dominant Celltype Accuracy') + xlab('Number of included RNA features sorted by variance')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying RNA features, effect on Cluster Accuracy')


outvecs$rna_features <- prettyNum(outvecs$rna_features, big.mark = ',')
outvecs$rna_features <- factor(outvecs$rna_features, levels = dfm$rna_features)


rna_acc_vln <- ggplot(outvecs, aes(rna_features, clusteraccuracy))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'steelblue', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  theme_light()+
  ylab('Cluster Dominant Celltype Accuracy') + xlab('Number of included RNA features sorted by variance')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying RNA features, effect on Cluster Accuracy')





pdf('outs/cerebralcortex/hvf_loop/SUMMARY_RNA_plot_ari_acc.pdf', height =3, width = 4)

rna_ari
rna_acc
rna_acc_vln

dev.off()




pdf('outs/cerebralcortex/hvf_loop/SUMMARY_RNA_table.pdf', height =8, width = 10)

SDAP::pdftable(dfs)

dev.off()





















### UPDATE, NOW REVERSE, HVF ATAC FEATURES



sobj <- readRDS("outs/cerebralcortex/sobj_binary_tfidf_lsi.rds")
ftt <- readRDS('outs/cerebralcortex/feature_types.rds')


### GET RNA FEATURES
## HVGs...
mat <- GetAssayData(sobj, assay = 'peaks', layer = 'counts')

#get RNA features
rna_ftt <- ftt[ftt$type == 'rna',]
mat <- mat[rownames(mat) %in% rna_ftt$feature,]


# find HVGs
vf <- FindVariableFeatures(mat)
vf <- vf[order(vf$vst.variance, decreasing = T),]


vf <- vf[1:2000,]
RNAVARFEATS <- rownames(vf)

rm(vf, mat, rna_ftt)
gc(full = T)


## HVGs...
mat <- GetAssayData(sobj, assay = 'peaks', layer = 'counts')

#get RNA features
table(ftt$type)
atac_ftt <- ftt[ftt$type == 'atac',]
mat <- mat[rownames(mat) %in% atac_ftt$feature,]


# # find HVGs
# vf <- FindVariableFeatures(mat)
# vf <- vf[order(vf$vst.variance, decreasing = T),]
# beepr::beep()

## ATAC: Find TOP features
library(Signac)
vf <- Signac::FindTopFeatures(mat)
vf <- vf[order(vf$percentile, decreasing = T),]



#add in amount of 
numhvgs = c(Inf, 95, 90, 80, 75, 60, 50, 25, 10, 0)
numhvgs = numhvgs/100


# save original, which used all features
sobj_orig <- sobj


nh = numhvgs[2]

timestart_atac = proc.time()





gc(full = T)
hvg_loop_res <- lapply(numhvgs, function(nh){
  
  
  
  hvf_filename <- paste0('outs/cerebralcortex/hvf_loop/', 'atac_', nh, '.rds')
  
  
  if(file.exists(hvf_filename)){
    
    
    message(nh, ' already done - reading in')
    
    return(readRDS(hvf_filename))
  } else{
    
    
    
    message('##################\n',
            '##################\n',
            '##################\n',
            '\n\n',
            'Now Running: ', nh, '\n',
            '\n\n',
            '##################\n',
            '##################\n',
            '##################\n'
            
    )
    
    
    
    #get sobj
    sobj <- sobj_orig
    
    gc(full = T)
    
    #get hvgs for this iteration
    vf_sub <- vf[vf$percentile > nh,,drop=F]
    
    #add in RNA VAR FEATURES, 2K for now
    varfeats <- c(rownames(vf_sub), RNAVARFEATS)
    
    VariableFeatures(sobj) <- varfeats
    
    
    
    #run tfidf and lsi
    mat <- GetAssayData(sobj, assay = 'peaks', layer = 'counts')
    mat <- mat[rownames(mat) %in% varfeats,]
    
    # featuere selection; or use all? 
    # sobj <- FindTopFeatures(sobj, min.cutoff = 'q0')
    # mat <- mat[rownames(mat) %in% VariableFeatures(sobj),]
    
    # tf idf
    tfidf <- Signac::RunTFIDF(mat)
    
    
    # lsi
    lsi <- RunSVD(tfidf)
    
    
    sobj[['lsi']] <- lsi
    
    
    rm(tfidf, lsi, mat)
    gc(full = T)
    
    
    ### downstream analysis: clustering and umap
    #umap
    sobj <- RunUMAP(sobj, reduction = 'lsi', dims=2:20, reduction.name = 'umap_lsi', reduction.key = 'umap_lsi')
    
    #graph
    sobj <- FindNeighbors(object = sobj, reduction = 'lsi', dims = 2:20)
    
    #clustering; try to get ~15 clusters
    # table(sobj$celltype)
    # length(table(sobj$celltype))
    sobj <- FindClusters(object = sobj, resolution = 0.2, graph.name = 'RNA_snn')
    
    
    
    # reduction = 'umap_lsi'
    # dp <- (DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'seurat_clusters') + SDAP::theme_dimplot()) / 
    #   (DimPlot(object = sobj, label = TRUE, repel = T, reduction = reduction, group.by = 'celltype') + SDAP::theme_dimplot())
    # 
    
    
    #get num atac and rna features
    ftt_in <- ftt[ftt$feature %in% varfeats,]
    ftt_in$type <- factor(ftt_in$type, levels = c('atac', 'rna'))
    ftt_tab <- table(ftt_in$type)
    
    # #prep title of patchwork
    # dptitle <- paste0( 'ATAC features: ', ftt_tab[1], '; ATAC quantile included: ', nh, '; RNA features: ', ftt_tab[2] )
    # 
    #try to get rand index, ari
    message('\n\nCALCULATE RAND INDEX\n\n\n')
    celltype_cluster <- sobj@meta.data[,c('celltype', 'seurat_clusters')]
    celltype_cluster$celltype <- factor(celltype_cluster$celltype)
    
    # rand <- fossil::rand.index(celltype_cluster[,1], celltype_cluster[,2])
    # adj.rand <- fossil::adj.rand.index(celltype_cluster[,1], celltype_cluster[,2])
    
    adj.rand <- round( pdfCluster::adj.rand.index(celltype_cluster[,1], celltype_cluster[,2]) , 2)
    
    
    
    
    ## also calculate "accuracy"
    twt <- as.matrix(table(celltype_cluster[,1], celltype_cluster[,2]))
    outs <- lapply(1:ncol(twt), function(i){
      x = twt[,i]
      
      (max(x) / sum(x)) * 100
    })
    outvec <- unlist(outs)
    names(outvec) <- colnames(twt)
    
    
    meanacc <- round(mean(outvec) , 2)
    minacc <- round(min(outvec) , 2)
    maxacc <- round( max(outvec) , 2)
    
    # dpsubtitle <- paste0('Adjusted Rand Index = ', adj.rand, '\n',
    #                      'Column Accuracy % [range] = ',  meanacc, ' [', minacc, ' - ', maxacc , ']')
    # 
    # 
    # 
    # dp <- dp + plot_annotation(title =  dptitle,
    #                            subtitle = dpsubtitle)
    
    
    #also get num clusts
    numclusts = length(levels(celltype_cluster[,2]))
    
    rm(mat, tfidf, lsi, vf_sub, ftt_in, varfeats, atac_features, celltype_cluster, rna_ftt, outs, twt,
       atac_ftt)
    gc(full = T)
    
    
    df = data.frame(atac_features = ftt_tab[1],
                    rna_features = ftt_tab[2],
                    ari = adj.rand,
                    meanacc = meanacc,
                    minacc = minacc,
                    maxacc = maxacc,
                    numclusts = numclusts)
    
    
    message('\nSaving...')
    
    saveRDS(list(sobj = sobj,
                 outvec = outvec,
                 df = df),
            file = hvf_filename
    )
    
    
    rm(sobj)
    gc(full = T)
    
    
    return(list(outvec = outvec,
                df = df)
    )
    
  } #file loop
  
})


gc(full = T)

beepr::beep()

hourspassed_atac <- (proc.time() - timestart_atac)[3]/60/60
names(hourspassed_atac) <- 'Hours'
hourspassed_atac

names(hvg_loop_res) <- paste0('atac_', numhvgs)


saveRDS(hvg_loop_res, "outs/cerebralcortex/hvf_loop/atac_sumstats.rds")
hvg_loop_res <- readRDS("outs/cerebralcortex/hvf_loop/atac_sumstats.rds")







dfs <- lapply(1:length(hvg_loop_res), function(i){
  hvg_loop_res[[i]]$df
})
dfs <- dplyr::bind_rows(dfs)
rownames(dfs) <- NULL

dfs$atac_percentile <- numhvgs
dfs <- dfs[,c(8,1:7)]

i = 2
outvecs = lapply(1:length(hvg_loop_res), function(i){
  thisres = hvg_loop_res[[i]]
  
  df = thisres$df
  ov = thisres$outvec
  
  quantile = names(hvg_loop_res)[i]
  quantile = gsub('atac_', '', quantile)
  
  if(quantile == 'Inf'){quantile = 'No ATAC features'}
  if(quantile == '0'){quantile = 'All ATAC features'}
  
  clusteraccdf = data.frame(atac_percentile = quantile,
                            cluster = names(ov),
                            clusteraccuracy = ov)
  
  clusteraccdf
  
  
})
outvecs <- dplyr::bind_rows(outvecs)
rownames(outvecs) <- NULL



## read in make some plots

dfs[1,1] <- 'No ATAC features'
dfs[nrow(dfs), 1] <- 'All ATAC features'

dfs_atac = dfs


write.csv(dfs_atac, 'outs/cerebralcortex/hvf_loop/ATAC_DF.csv', row.names = F, quote = F)

dfm <- dfs_atac
dfm$atac_percentile <- factor(dfm$atac_percentile , levels = dfm$atac_percentile )



atac_ari <- ggplot(dfm, aes(atac_percentile, ari))+
  geom_point()+
  theme_light()+
  ylab('Adjusted Rand Index') + xlab('Quantile of top ATAC features included')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying ATAC features, effect on A.R.I.')


atac_acc <- ggplot(dfm, aes(atac_percentile, meanacc))+
  geom_point()+
  theme_light()+
  ylab('Mean Cluster Dominant Celltype Accuracy') + xlab('Quantile of top ATAC features included')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying ATAC features, effect on Cluster Accuracy')


outvecs$atac_percentile <- factor(outvecs$atac_percentile, levels = dfm$atac_percentile)

atac_acc_vln <- ggplot(outvecs, aes(atac_percentile, clusteraccuracy))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'firebrick', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  theme_light()+
  ylab('Cluster Dominant Celltype Accuracy') + xlab('Quantile of top ATAC features included')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying ATAC features, effect on Cluster Accuracy')






pdf('outs/cerebralcortex/hvf_loop/SUMMARY_ATAC_plot_ari_acc.pdf', height =3, width = 4)

atac_ari
atac_acc
atac_acc_vln

dev.off()




pdf('outs/cerebralcortex/hvf_loop/SUMMARY_ATAC_table.pdf', height =8, width = 10)

SDAP::pdftable(dfs_atac)

dev.off()











#### read in sweep, make a umap of first, mid, last, to highlight some celltypes


#for RNA
rna0 <- 'outs/cerebralcortex/hvf_loop/rna_0.rds'
rna15000 <- 'outs/cerebralcortex/hvf_loop/rna_15000.rds'
rnaall <- 'outs/cerebralcortex/hvf_loop/rna_30201.rds'

plist = list(rna0=rna0, rna15000=rna15000,rnaall=rnaall)




labsize = 4


dp_l <- lapply(plist, function(inpath){
  
  message(inpath)
  inlist <- readRDS(inpath)
  
  # inlist <- rna0
  sobj <- inlist$sobj
  df <- inlist$df
  
  d1 <- DimPlot(sobj, group.by = 'RNA_snn_res.0.2', reduction = 'umap_lsi', label = T, repel = T, label.size = labsize) + SDAP::theme_dimplot()
  d2 <- DimPlot(sobj, group.by = 'celltype', reduction = 'umap_lsi',  label = T, repel = 2, label.size = labsize)+ SDAP::theme_dimplot()
  
  
  dptitle <- paste0( 'ATAC features: ', df[1,1], '; RNA features: ', df[1,2])
  dpsubtitle <- paste0('Adjusted Rand Index = ', df[1,3], '\n',
                       'Column Accuracy % [range] = ',  df[1,4], ' [', df[1,5], ' - ', df[1,6] , ']')
  
  
  
  d = d1/d2
  d <- d + plot_annotation(title =  dptitle,
                           subtitle = dpsubtitle)
  
  
  return(d)
  
  
})

# dir.create('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/')
# pdf('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/RNA.pdf', height = 9, width = 7)
# 
# dp_l
# dev.off()


#print one at a time...
lapply(1:length(dp_l), function(i){
  
  
  dp = dp_l[[i]]
  dpname = names(dp_l)[i]
  
  pdfname = paste0('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/RNA_', dpname, '.pdf')
  pdf(pdfname, height = 9, width = 7)
  
  
  print(dp)
  
  dev.off()
  
  return(dpname)
  
})




### same for atac

atac_none <- 'outs/cerebralcortex/hvf_loop/atac_Inf.rds'
atac_0.95 <- 'outs/cerebralcortex/hvf_loop/atac_0.95.rds'
atac0.5 <- 'outs/cerebralcortex/hvf_loop/atac_0.5.rds'
atac_all <- 'outs/cerebralcortex/hvf_loop/atac_0.rds'


plist = list(atac_none=atac_none,atac_0.95= atac_0.95, atac0.5=atac0.5,atac_all=atac_all)
inpath = plist[[3]]





dp_l <- lapply(plist, function(inpath){
  
  message(inpath)
  inlist <- readRDS(inpath)
  
  # inlist <- rna0
  sobj <- inlist$sobj
  df <- inlist$df
  
  d1 <- DimPlot(sobj, group.by = 'RNA_snn_res.0.2', reduction = 'umap_lsi', label = T, repel = T, label.size = labsize) + SDAP::theme_dimplot()
  d2 <- DimPlot(sobj, group.by = 'celltype', reduction = 'umap_lsi',  label = T, repel = 2, label.size = labsize)+ SDAP::theme_dimplot()
  
  
  dptitle <- paste0( 'ATAC features: ', df[1,1], '; RNA features: ', df[1,2])
  dpsubtitle <- paste0('Adjusted Rand Index = ', df[1,3], '\n',
                       'Column Accuracy % [range] = ',  df[1,4], ' [', df[1,5], ' - ', df[1,6] , ']')
  
  
  
  d = d1/d2
  d <- d + plot_annotation(title =  dptitle,
                           subtitle = dpsubtitle)
  
  
  return(d)
  
  
})


# dir.create('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/')
# pdf('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/ATAC.pdf', height = 9, width = 7)
# 
# dp_l
# dev.off()




#print one at a time...
lapply(1:length(dp_l), function(i){
  
  
  dp = dp_l[[i]]
  dpname = names(dp_l)[i]
  
  pdfname = paste0('outs/cerebralcortex/hvf_loop/S10_highlight_improvements/ATAC_', dpname, '.pdf')
  pdf(pdfname, height = 9, width = 7)
  
  
  print(dp)
  
  dev.off()
  
  return(dpname)
  
})





