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




#read object and feature types (rna vs atac)
sobj <- readRDS('outs/cerebralcortex/sobj_binary_tfidf_lsi.rds')
ftt <- readRDS('outs/cerebralcortex/feature_types.rds')


### FIGURES


# Dimplot of UMAP with Clusters
# Dimplot of UMAP with published celltypes
# Alluvial Plot of Celltypes -> Clusters
# Heatmap Matrix of Celltypes to clusters



#set heights/ widths
dimplot_h <- 4
dimplot_w <- dimplot_h

dotplot_h <- 5
dotplot_w <- 8

matrix_h <- 6
matrix_w <- 7





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










##########################################################################################
##########################################################################################
##########################################################################################

# Dimplot of UMAP with Clusters

##########################################################################################
##########################################################################################
##########################################################################################




pa <- DimPlot(sobj, group.by = "RNA_snn_res.0.2", reduction = 'umap_lsi', 
              label = T, repel = T, label.size = 2) + 
  ggtitle('Binarized with TF/IDF and LSI, Louvain resolution = 0.2')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot()+
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))


pdf('outs/cerebralcortex//figures/panel_A.pdf', height = dimplot_h, width = dimplot_w)
pa
pa + ggtitle('')
pa +NoLegend()
dev.off()


svg('outs/cerebralcortex/figures/panel_A.svg', height = dimplot_h, width = dimplot_w)
pa +NoLegend()
dev.off()









##########################################################################################
##########################################################################################
##########################################################################################

# Dimplot of UMAP with published celltypes

##########################################################################################
##########################################################################################
##########################################################################################




pb <- DimPlot(sobj, group.by = "celltype", reduction = 'umap_lsi', 
              label = T, repel = T, label.size = 2) + 
  ggtitle('Binarized with TF/IDF and LSI, Published Celltypes')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot()+
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))


pdf('outs/cerebralcortex//figures/panel_B.pdf', height = dimplot_h, width = dimplot_w)
pb
pb + ggtitle('')
pb +NoLegend()
dev.off()


svg('outs/cerebralcortex/figures/panel_B.svg', height = dimplot_h, width = dimplot_w)
pb +NoLegend()
dev.off()



#remove cell label in UMAP
pb <- DimPlot(sobj, group.by = "celltype", reduction = 'umap_lsi', 
              label = F, repel = T, label.size = 2) + 
  ggtitle('Binarized with TF/IDF and LSI, Published Celltypes')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot()+
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))


pdf('outs/cerebralcortex//figures/panel_B_nocelllab.pdf', height = dimplot_h, width = dimplot_w)
pb
pb + ggtitle('')
pb +NoLegend()
dev.off()


svg('outs/cerebralcortex/figures/panel_B_nocelllab.svg', height = dimplot_h, width = dimplot_w)
pb +NoLegend()
dev.off()






##########################################################################################
##########################################################################################
##########################################################################################

# Alluvial Plot of Celltypes -> Clusters

##########################################################################################
##########################################################################################
##########################################################################################


labelsdf <- sobj@meta.data[,c('celltype', 'RNA_snn_res.0.2')]
colnames(labelsdf)[1] = 'Published Celltype'


pc <- SDAP::alluvialplot(labelsdf = labelsdf, repel = T, size = 1.5)

dotplot_h = 5
dotplot_w = 6


pdf('outs/cerebralcortex//figures/panel_C.pdf', height = dotplot_h, width = dotplot_w)
pc
pc + theme(legend.position = 'bottom')
pc + NoLegend()
dev.off()


pc <- SDAP::alluvialplot(labelsdf = labelsdf, repel = T, size = 5)
pdf('outs/cerebralcortex//figures/panel_C_legendbottom.pdf', height = 10, width = 10)

pc + theme(legend.position = 'bottom')

dev.off()


svg('outs/cerebralcortex/figures/panel_C.svg', height = dimplot_h, width = dimplot_w)
pc
dev.off()









##########################################################################################
##########################################################################################
##########################################################################################

# Heatmap Matrix of Celltypes to clusters

##########################################################################################
##########################################################################################
##########################################################################################




### map clusters to top celltype...
celltype_cluster <- sobj@meta.data[,c('celltype', 'seurat_clusters')]
cmax <- get_top_celltype(celltype_cluster)

#sort cmax alphabetically while retaining the vecotr names which are the cluster numbers
cmax <- cmax[str_order(cmax, numeric = T)]

#using the sorted cluster mapped labels, reorder factor levels...
sobj$clusters_sorted_by_top_celltype <- sobj$seurat_clusters
sobj$clusters_sorted_by_top_celltype <- factor(sobj$clusters_sorted_by_top_celltype, levels = names(cmax))


#make a two way table and heatmap
twt <- as.matrix(table(sobj$celltype, sobj$clusters_sorted_by_top_celltype))

twt_scale <- scale(twt)

library(ComplexHeatmap)
library(grid)
pd <- ComplexHeatmap::Heatmap(twt_scale, cluster_rows = F, cluster_columns = F, 
                              rect_gp = grid::gpar(col = "white", lwd = 0.5),
                              border_gp = grid::gpar(col = "black", lwd = 2),
                              column_title = 'Binary TF/IDF LSI-derived Clusters', column_title_side = 'bottom',
                              column_names_rot = 45,
                              column_names_gp = grid::gpar(fontsize = 12),
                              row_names_gp = grid::gpar(fontsize = 8),
                              row_title = 'Published Celltypes',
                              row_names_side = 'left',
                              show_heatmap_legend = F,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid::grid.text(sprintf("%.0f", twt[i, j]), x, y, gp = gpar(fontsize = 7, col = 'white'))
                              })


pdf('outs/cerebralcortex//figures/panel_D.pdf', height = matrix_h, width = matrix_w)
pd
dev.off()


svg('outs/cerebralcortex/figures/panel_D.svg', height = matrix_h, width = matrix_w)
pd
dev.off()









# check accuracy

# binary tf/idf lsi

twt <- as.matrix(table(sobj$celltype, sobj$clusters_sorted_by_top_celltype))
outs <- lapply(1:ncol(twt), function(i){
  x = twt[,i]
  
  (max(x) / sum(x)) * 100
})
outvec <- unlist(outs)
names(outvec) <- colnames(twt)
outvec

sort(outvec)


mean(outvec)
min(outvec) 
max(outvec)

# > mean(outvec)
# [1] 86.93417
# > min(outvec) 
# [1] 22.95174
# > max(outvec)
# [1] 99.21671

accvec <- outvec
sizes <- colSums(twt)

df <- data.frame(accvec,sizes, cluster=names(accvec))
df <- df[order(df$sizes),]


ggplot(df, aes(sizes,accvec, label=cluster))+
  geom_point()+
  ggrepel::geom_text_repel(size = 5)+
  xlab('cluster size') + ylab('accuracy')








###### cluster 10 mixed cluster markers ######

sobj <- readRDS("outs/cerebralcortex/sobj_nonbinary_rnaonly.rds")



#use makrers from fig 1F of paper
markers <- c('VIM', 'PAX6', 'HES5', 'EOMES', 'SATB2', 'SLC17A7', 'NEUROD2', 'GAD2', 'LHX6', 'VIP', 'OLIG1', 'SOX10', 'AQP4', 'OPALIN', 'PTPRC', 'CLDN5', 'PDGFRB', 'COL1A2')



sobj$celltype_dotplot <- factor(sobj$celltype, levels = c('RG','IPC', 'EN-fetal-early', 'EN-fetal-late', 'EN', 'IN-fetal', 'IN-MGE', 'IN-CGE', 'OPC', 'Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Pericytes', 'VSMC'))
dp_ctmarkers_likepaper <- DotPlot(sobj, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 

celltype_cluster <- sobj@meta.data[,c('celltype_dotplot', 'binary_RNA_snn_res.0.2')]
cmax = get_top_celltype(celltype_cluster)

#try to match levs
levs <- levels(sobj$celltype_dotplot)
levs <- levs[levs%in%cmax]
cmax = factor(cmax, levels = levs)
cmax = cmax[order(cmax)]

sobj$binary_RNA_snn_res.0.2_cmax <- factor(sobj$binary_RNA_snn_res.0.2, levels = names(cmax))

dp_ctmarkers_binclusters <- DotPlot(sobj, features = markers, group.by = 'binary_RNA_snn_res.0.2_cmax')+
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

vlnlot_qc <- VlnPlot(sobj, c('nCount_RNA', 'nFeature_RNA', 'percent_mt'), ncol = 1, pt.size = 0.1) + NoLegend()



## make two more dot plots
# 1. with all cells except cluster 10, using the paper markers, and using the paper cell types
md <- sobj@meta.data
md <- md[md$binary_RNA_snn_res.0.2 != 10,]
sobj_noc10 <- sobj[,rownames(md)]


dp_noc10 <- DotPlot(sobj_noc10, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) +
  ggtitle('No Cluster 10 cells')


# 2. with only cluster 10 cells, using paper markers, using paper cell types
md <- sobj@meta.data
md <- md[md$binary_RNA_snn_res.0.2 == 10,]
sobj_onlyc10 <- sobj[,rownames(md)]


dp_onlyc10 <- DotPlot(sobj_onlyc10, features = markers, group.by = 'celltype_dotplot')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) +
  ggtitle('Only Cluster 10 cells')




## plot top 5 markers per cluster
m <- read.csv("outs/cerebralcortex/markers_binaryclusters_RNA_snn_res.0.2_calc-from-nonbinary.csv")


n = 3

m$score <- (m$pct.1 - m$pct.2) * m$avg_log2FC
m <- m[order(m$score, decreasing = T),] #sort by lfc and remove; this way top marker is unique cluster-wise
m <- m[!duplicated(m$gene),]
m <- m[order(m$cluster),]
top <- m %>% group_by(cluster) %>% top_n(n = n, wt = score)
genes <- unique(top$gene)

sobj <- GetResidual(sobj, features = genes)


dp_binclust_markers <- DotPlot(sobj, features = rev(genes), group.by = 'binary_RNA_snn_res.0.2')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 5),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) +
  ggtitle('Binary Cluster markers (computed from non-binarized data)')





## dotplot like paper
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_dotplot_likepaper"
pdf( paste0(plotname, '.pdf'), 
    height =dotplot_h, width = dotplot_w)

dp_ctmarkers_likepaper

dev.off()

svg( paste0(plotname, '.svg') , 
    height =dotplot_h, width = dotplot_w)

dp_ctmarkers_likepaper
dp_ctmarkers_likepaper+NoLegend()

dev.off()


## dotplot of paper markers in our clusters
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_dotplot_BinaryClustersSortedByDominantCelltype_nonBinarizedData"
pdf( paste0(plotname, '.pdf'), 
     height =dotplot_h, width = dotplot_w)

dp_ctmarkers_binclusters

dev.off()

svg( paste0(plotname, '.svg') , 
     height =dotplot_h, width = dotplot_w)

dp_ctmarkers_binclusters
dp_ctmarkers_binclusters+NoLegend()

dev.off()



## vln plot of qc
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_vln_qc"
pdf( paste0(plotname, '.pdf'), 
     height =dotplot_h, width = dotplot_w/1.5)

vlnlot_qc

dev.off()

svg( paste0(plotname, '.svg') , 
     height =dotplot_h, width = dotplot_w/1.5)

vlnlot_qc

dev.off()


## dot plot w/o c10
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_dotplot_noC10Cells"
pdf( paste0(plotname, '.pdf'), 
     height =dotplot_h, width = dotplot_w)

dp_noc10
dp_noc10+NoLegend()

dev.off()

svg( paste0(plotname, '.svg') , 
     height =dotplot_h, width = dotplot_w)

dp_noc10

dev.off()




## dot plot w only c10
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_dotplot_onlyC10Cells"
pdf( paste0(plotname, '.pdf'), 
     height =dotplot_h, width = dotplot_w)

dp_onlyc10
dp_onlyc10+NoLegend()

dev.off()

svg( paste0(plotname, '.svg') , 
     height =dotplot_h, width = dotplot_w)

dp_onlyc10

dev.off()



## dotplot plot w new markers
plotname = "outs/cerebralcortex/figures/C10IdentityPlots_dotplot_BinaryMarkers_ComputedFromNonBinary"
pdf( paste0(plotname, '.pdf'), 
     height =dotplot_h, width = dotplot_w)

dp_binclust_markers
dp_binclust_markers+NoLegend()

dev.off()

svg( paste0(plotname, '.svg') , 
     height =dotplot_h, width = dotplot_w)

dp_binclust_markers

dev.off()
