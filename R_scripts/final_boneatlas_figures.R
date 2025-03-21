library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(ComplexHeatmap)
library(grid)


set.seed(2022)




### BINARIZATION OF DATA FROM BARYAWNO ET AL 2019
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6570562/

setwd("~/Dropbox/Result_from_Alex/deyoudata/binarizationproject/")

sobjfull <- readRDS("outs/boneatlas/sobjfull.rds")
sobj <- readRDS('outs/boneatlas/sobj_binary.rds')
sobjlsi <- readRDS('outs/boneatlas/sobj_tfidf_lsi.rds')

pubm <- readRDS('outs/boneatlas/PARSED_MARKERS_NIHMS1529101-supplement-8.rds')


cellmaplist <- readRDS('outs/cellmaplist.rds')










### FIGURES


# Panel A: umap colored by celltype (non-binary) 
# Panel B: Bubble plot of their markers in Celltypes (maybe panel B) (try to pick the paper's canonical markers)
# Panel C heatmap/table of Celltype vs TF/IDF clusters
# Panel D: Umap colored by binary clusters
# Panel E: Bubble plot of their makrers in binarized clusters same markers as from B)
# Panel F heatmap/table of Celltype vs clusters
# Panel G: Umap colored by binary TF/IDF clusters
# Panel H: Bubble plot of their makrers in binary TF/IDF clusters (same markers as from B)
# Panel I heatmap/table of Celltype vs TF/IDF clusters



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

# A. UMAP of the non-binarized data, colored by CLUSTERS.

##########################################################################################
##########################################################################################
##########################################################################################



#sort levels of celltypes
# ctlevs <-  c('Lepr-MSC', 'OLC-1', 'OLC-2', 'Fibro-1', 'Fibro-2', 'Fibro-3', 'Fibro-4', 'Fibro-5', 'Chondro', 'Chondro-hyper', 'Chondro-prehyper-2', 'Chondro-progen', 'Chondro-prol/rest', 'Pericytes', 'EC-arterial', 'EC-arteriolar', 'EC-Sinusoidal')
# sobjfull$Celltype <- factor(sobjfull$Celltype, levels = ctlevs)

pa <- DimPlot(sobjfull, group.by = 'cluster_celltype', reduction = 'umap', 
              label = T, repel = T, label.size = 2) + 
  ggtitle('Non-Binarized Data, Louvain resolution 0.5')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot()+
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))


pdf('outs/boneatlas/figures/panel_A.pdf', height = dimplot_h, width = dimplot_w)
pa
pa +NoLegend()
dev.off()


svg('outs/boneatlas/figures/panel_A.svg', height = dimplot_h, width = dimplot_w)
pa +NoLegend()
dev.off()








##########################################################################################
##########################################################################################
##########################################################################################

# B. Bubble plot showing the expression of the cell type markers across clusters in A. 

##########################################################################################
##########################################################################################
##########################################################################################



#plot top 2 markers

n = 2
# top <- pubm %>% group_by(celltype) %>% top_n(n = n, wt = avg_logFC)
# get celltypes in...
topcts <- get_top_celltype(sobjfull@meta.data[,c('Celltype', 'seurat_clusters')])

pubm <- pubm[pubm$celltype %in% topcts,]



# just get all genes, but put them in alphabetical cell type order
pubm$avg_logFC <- as.numeric(pubm$avg_logFC)
pubm <- pubm[pubm$avg_logFC > 0,]
pubm_all <- pubm
cts <- levels(pubm_all$celltype)
cts <- str_sort(cts, numeric = T)
pubm_all$celltype <- factor(pubm_all$celltype, levels = cts)

#make the scoring thing..
pubm_all$score <- (pubm_all$pct.1 - pubm_all$pct.2) * pubm_all$avg_logFC
pubm_all <- pubm_all[order(pubm_all$score, decreasing = T),] #sort by lfc and remove; this way top marker is unique cluster-wise
pubm_all <- pubm_all[!duplicated(pubm_all$gene),]
pubm_all <- pubm_all[order(pubm_all$celltype),]

n = 2
top <- pubm_all %>% group_by(celltype) %>% top_n(n = n, wt = score)
genes <- unique(top$gene)

sobjfull <- GetResidual(sobjfull, features = genes)






sobjfull$cluster_celltype_cmaxsort_rev <- factor(sobjfull$cluster_celltype_cmaxsort, levels = rev(levels(sobjfull$cluster_celltype_cmaxsort)))

# change "cluster" to "C"
levs <- levels(sobjfull$cluster_celltype_cmaxsort_rev )
nulevs <- gsub('Cluster_', 'C', levs)
sobjfull$cluster_celltype_cmaxsort_rev <- plyr::mapvalues(sobjfull$cluster_celltype_cmaxsort_rev,
                                                          levs, nulevs)



pb <- DotPlot(sobjfull, features = genes, group.by = 'cluster_celltype_cmaxsort_rev', assay = 'SCT') + 
  #coord_flip() + 
  xlab('') + ylab('')+
  # ggtitle('Non-Binarized Original Dataset', subtitle = 'Markers from Publication')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 



pdf('outs/boneatlas/figures/panel_B.pdf', height = dotplot_h, width = dotplot_w)
pb
pb+NoLegend()
dev.off()


svg('outs/boneatlas/figures/panel_B.svg', height = dotplot_h, width = dotplot_w)
pb
dev.off()









##########################################################################################
##########################################################################################
##########################################################################################

# C -Panel C heatmap/table of Celltype vs full clusters


##########################################################################################
##########################################################################################
##########################################################################################


sobjfull$Celltype_abc <- factor(sobjfull$Celltype, levels = str_sort(levels(sobjfull$Celltype), numeric = T) )

#for cluster_celltype_cmaxsort - replace "Cluster_" with C"
# change "cluster" to "C"
levs <- levels(sobjfull$cluster_celltype_cmaxsort )
nulevs <- gsub('Cluster_', 'C', levs)
sobjfull$cluster_celltype_cmaxsort_C <- plyr::mapvalues(sobjfull$cluster_celltype_cmaxsort,
                                                          levs, nulevs)

#make a two way table and heatmap

twt <- as.matrix(table(sobjfull$Celltype_abc, sobjfull$cluster_celltype_cmaxsort_C))

twt_scale <- scale(twt)

library(ComplexHeatmap)
library(grid)
pc <- ComplexHeatmap::Heatmap(twt_scale, cluster_rows = F, cluster_columns = F, 
                              rect_gp = grid::gpar(col = "white", lwd = 0.5),
                              border_gp = grid::gpar(col = "black", lwd = 2),
                              column_title = 'Non-Binarized Clusters', column_title_side = 'bottom',
                              column_names_rot = 45,
                              column_names_gp = grid::gpar(fontsize = 8),
                              row_names_gp = grid::gpar(fontsize = 8),
                              row_title = "Published Celltypes",
                              row_names_side = 'left',
                              show_heatmap_legend = F,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid::grid.text(sprintf("%.0f", twt[i, j]), x, y, gp = gpar(fontsize = 7, col = 'white'))
                              })



pdf('outs/boneatlas/figures/panel_C.pdf', height = matrix_h, width = matrix_w)
pc
dev.off()


svg('outs/boneatlas/figures/panel_C.svg', height = matrix_h, width = matrix_w)
pc
dev.off()
















##########################################################################################
##########################################################################################
##########################################################################################

# # Panel D: Umap colored by binary clusters

##########################################################################################
##########################################################################################
##########################################################################################


pd <- DimPlot(sobj, group.by = 'seurat_clusters', reduction = 'umap', 
              label = T, repel = T, label.size = 2) + 
  ggtitle('Binarized, Louvain resolution 0.5')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot()+ 
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))


pdf('outs/boneatlas/figures/panel_D.pdf', height = dimplot_h, width = dimplot_w)
pd
pd+
  NoLegend()
dev.off()


svg('outs/boneatlas/figures/panel_D.svg', height = dimplot_h, width = dimplot_w)
pd+
  NoLegend()
dev.off()



##########################################################################################
##########################################################################################
##########################################################################################

# Panel E: Bubble plot of their makrers in binarized clusters same markers as from B)


##########################################################################################
##########################################################################################
##########################################################################################


# # get first celtype from each cluster
# pubm_bin <- pubm
# cts <- levels(sobj$cluster_celltype_ctsort)
# cts <- str_split_fixed(cts, '--', n=2)[,2]
# cts <- str_split_fixed(cts, '_', n=2)[,1]
# pubm_bin <- pubm_bin[pubm_bin$celltype %in% cts,]
# pubm_bin$celltype <- factor(pubm_bin$celltype, levels = unique(cts))
# pubm_bin <- pubm_bin[!duplicated(pubm_bin$gene),]
# pubm_bin <- pubm_bin[order(pubm_bin$celltype),]
# 
# top <- pubm_bin %>% group_by(celltype) %>% top_n(n = n, wt = avg_logFC)
# 
# genes_bin <- unique(top$gene)

#use same genes


sobjfull$binary_cluster_celltype_cmaxsort_rev <- factor(sobj$cluster_celltype_cmaxsort, levels = rev(levels(sobj$cluster_celltype_cmaxsort) ))
levs <- levels(sobjfull$binary_cluster_celltype_cmaxsort_rev)
levs <- gsub('Cluster_', '', levs)
levs <- str_split_fixed(levs, '--', 2)[,1]
sobjfull$binary_cluster_celltype_cmaxsort_rev <- plyr::mapvalues(sobjfull$binary_cluster_celltype_cmaxsort_rev ,
                                                                 from = levels(sobjfull$binary_cluster_celltype_cmaxsort_rev ),
                                                                 to = levs)


pe <- DotPlot(sobjfull, features = genes, group.by = 'binary_cluster_celltype_cmaxsort_rev', assay = 'SCT') + 
  xlab('') + ylab('')+
  # coord_flip()+
  # ggtitle('DotPlot of Markers in Binarized Data', subtitle = 'Markers from Publication')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        # axis.text.y = element_text(size = 6),
        legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 


pdf('outs/boneatlas/figures/panel_E.pdf', height = dotplot_h, width = dotplot_w)
pe
pe+NoLegend()
dev.off()


svg('outs/boneatlas/figures/panel_E.svg', height = dotplot_h, width = dotplot_w)
pe
dev.off()







##########################################################################################
##########################################################################################
##########################################################################################

# Panel F heatmap/table of full clusters vs binary clusters


##########################################################################################
##########################################################################################
##########################################################################################


sobjfull$binary_cluster_celltype_cmaxsort <- sobj$cluster_celltype_cmaxsort

levs <- levels(sobjfull$binary_cluster_celltype_cmaxsort)
levs <- gsub('Cluster_', '', levs)
levs <- str_split_fixed(levs, '--', 2)[,1]
sobjfull$binary_cluster_celltype_cmaxsort <- plyr::mapvalues(sobjfull$binary_cluster_celltype_cmaxsort ,
                                                             from = levels(sobjfull$binary_cluster_celltype_cmaxsort ),
                                                             to = levs)


#make a two way table and heatmap
twt <- as.matrix(table(sobjfull$cluster_celltype_cmaxsort_C, sobjfull$binary_cluster_celltype_cmaxsort))

twt_scale <- scale(twt)

library(ComplexHeatmap)
library(grid)
pf <- ComplexHeatmap::Heatmap(twt_scale, cluster_rows = F, cluster_columns = F, 
                              rect_gp = grid::gpar(col = "white", lwd = 0.5),
                              border_gp = grid::gpar(col = "black", lwd = 2),
                              column_title = 'Binary-derived Clusters', column_title_side = 'bottom',
                              column_names_rot = 45,
                              column_names_gp = grid::gpar(fontsize = 12),
                              row_names_gp = grid::gpar(fontsize = 8),
                              row_title = 'Non-Binarized Clusters',
                              row_names_side = 'left',
                              show_heatmap_legend = F,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid::grid.text(sprintf("%.0f", twt[i, j]), x, y, gp = gpar(fontsize = 7, col = 'white'))
                              })



pdf('outs/boneatlas/figures/panel_F.pdf', height = matrix_h, width = matrix_w)
pf
dev.off()


svg('outs/boneatlas/figures/panel_F.svg', height = matrix_h, width = matrix_w)
pf
dev.off()









##########################################################################################
##########################################################################################
##########################################################################################

# panel G UMAP of TFIDF


##########################################################################################
##########################################################################################
##########################################################################################


# ctlevs <-  c('Lepr-MSC', 'OLC-1', 'OLC-2', 'Fibro-1', 'Fibro-2', 'Fibro-3', 'Fibro-4', 'Fibro-5', 'Chondro', 'Chondro-hyper', 'Chondro-prehyper-2', 'Chondro-progen', 'Chondro-prol/rest', 'Pericytes', 'EC-arterial', 'EC-arteriolar', 'EC-Sinusoidal')
# sobjlsi$Celltype <- factor(sobj$Celltype, levels = ctlevs)

pg <- DimPlot(sobjlsi, reduction = 'umap_lsi', 
              group.by = 'seurat_clusters',label = T, repel = T,
              label.size = 2) + 
  ggtitle('Binarized with TF/IDF and LSI, Louvain resolution 0.5')+
  xlab('UMAP 1') + ylab('UMAP 2')+ 
  SDAP::theme_dimplot() + 
  theme(plot.title = element_text(hjust = 0, 
                                  vjust = 1, face = "bold", size = 10))



pdf('outs/boneatlas/figures/panel_G.pdf', height = dimplot_h, width = dimplot_w)
pg
pg+
  NoLegend()

dev.off()


svg('outs/boneatlas/figures/panel_G.svg', height = dimplot_h, width = dimplot_w)
pg+
  NoLegend()
dev.off()





##########################################################################################
##########################################################################################
##########################################################################################

# H - Bubble plot showing the expression of cell type markers across clusters in E.


##########################################################################################
##########################################################################################
##########################################################################################



# # get first celtype from each cluster
# pubm_lsi <- pubm
# cts <- levels(sobjlsi$cluster_celltype_ctsort)
# cts <- str_split_fixed(cts, '--', n=2)[,2]
# cts <- str_split_fixed(cts, '_', n=2)[,1]
# pubm_lsi <- pubm_lsi[pubm_lsi$celltype %in% cts,]
# pubm_lsi$celltype <- factor(pubm_lsi$celltype, levels = unique(cts))
# pubm_lsi <- pubm_lsi[!duplicated(pubm_lsi$gene),]
# pubm_lsi <- pubm_lsi[order(pubm_lsi$celltype),]
# 
# top <- pubm_lsi %>% group_by(celltype) %>% top_n(n = n, wt = avg_logFC)
# 
# genes_lsi <- unique(top$gene)
# sobjfull <- GetResidual(sobjfull, features = genes_lsi)


sobjfull$cluster_celltype_cmaxsort_rev <- factor(sobjlsi$cluster_celltype_cmaxsort, levels = rev(levels(sobjlsi$cluster_celltype_cmaxsort) ))
levs <- levels(sobjfull$cluster_celltype_cmaxsort_rev)
levs <- gsub('Cluster_', '', levs)
levs <- str_split_fixed(levs, '--', 2)[,1]
sobjfull$cluster_celltype_cmaxsort_rev <- plyr::mapvalues(sobjfull$cluster_celltype_cmaxsort_rev ,
                                                          from = levels(sobjfull$cluster_celltype_cmaxsort_rev ),
                                                          to = levs)




# sobjfull$cluster_celltype_cmaxsort_rev <- factor(sobjlsi$cluster_celltype_cmaxsort, levels = rev(levels(sobjlsi$cluster_celltype_cmaxsort) ))


ph <- DotPlot(sobjfull, features = genes, group.by = 'cluster_celltype_cmaxsort_rev', assay = 'SCT') + 
  xlab('') + ylab('')+
  # coord_flip()+
  # ggtitle('DotPlot of Markers in Binarized Data', subtitle = 'Markers from Publication')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        # axis.text.y = element_text(size = 6),
        # legend.title=element_text(size=7),
        plot.margin = unit(c(0,0,0,0),"cm")
  ) 




pdf('outs/boneatlas/figures/panel_H.pdf', height = dotplot_h, width = dotplot_w)
ph
ph + NoLegend()
dev.off()


svg('outs/boneatlas/figures/panel_H.svg', height = dotplot_h, width = dotplot_w)
ph
dev.off()
























##########################################################################################
##########################################################################################
##########################################################################################

# I - heatmap of LSI versus Full clusters


##########################################################################################
##########################################################################################
##########################################################################################



# heatmap
sobjfull$lsi_cluster_celltype_cmaxsort <- sobjlsi$cluster_celltype_cmaxsort



levs <- levels(sobjfull$lsi_cluster_celltype_cmaxsort)
levs <- gsub('Cluster_', '', levs)
levs <- str_split_fixed(levs, '--', 2)[,1]
sobjfull$lsi_cluster_celltype_cmaxsort <- plyr::mapvalues(sobjfull$lsi_cluster_celltype_cmaxsort ,
                                                          from = levels(sobjfull$lsi_cluster_celltype_cmaxsort ),
                                                          to = levs)



#make a two way table and heatmap
twt <- as.matrix(table(sobjfull$cluster_celltype_cmaxsort_C, sobjfull$lsi_cluster_celltype_cmaxsort))

twt_scale <- scale(twt)

library(ComplexHeatmap)
library(grid)
pi <- ComplexHeatmap::Heatmap(twt_scale, cluster_rows = F, cluster_columns = F, 
                              rect_gp = grid::gpar(col = "white", lwd = 0.5),
                              border_gp = grid::gpar(col = "black", lwd = 2),
                              column_title = 'Binary TF/IDF LSI-derived Clusters', column_title_side = 'bottom',
                              column_names_rot = 45,
                              column_names_gp = grid::gpar(fontsize = 12),
                              row_names_gp = grid::gpar(fontsize = 8),
                              row_title = 'Non-Binarized Clusters',
                              row_names_side = 'left',
                              show_heatmap_legend = F,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid::grid.text(sprintf("%.0f", twt[i, j]), x, y, gp = gpar(fontsize = 7, col = 'white'))
                              })


pdf('outs/boneatlas/figures/panel_I.pdf', height = matrix_h, width = matrix_w)
pi
dev.off()


svg('outs/boneatlas/figures/panel_I.svg', height = matrix_h, width = matrix_w)
pi
dev.off()






### MEAN ACCURACY ###


## full clusters vs celltype

#binary
twt <- as.matrix(table(sobjfull$Celltype, sobjfull$cluster_celltype))
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

accvec <- outvec
sizes <- colSums(twt)

df <- data.frame(accvec,sizes, cluster=names(accvec))
df <- df[order(df$sizes),]


ggplot(df, aes(sizes,accvec, label=cluster))+
  geom_point()+
  ggrepel::geom_text_repel(size = 1.5)+
  xlab('cluster size') + ylab('accuracy')



#binary
twt <- as.matrix(table(sobjfull$cluster_celltype_ctsort, sobjfull$binary_cluster_celltype_cmaxsort))
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

accvec <- outvec
sizes <- colSums(twt)

df <- data.frame(accvec,sizes, cluster=names(accvec))
df <- df[order(df$sizes),]


ggplot(df, aes(sizes,accvec, label=cluster))+
  geom_point()+
  ggrepel::geom_text_repel(size = 1.5)+
  xlab('cluster size') + ylab('accuracy')




# tf/idf and lsi
twt <- as.matrix(table(sobjfull$cluster_celltype_ctsort, sobjfull$lsi_cluster_celltype_cmaxsort))

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

accvec <- outvec
sizes <- colSums(twt)

df <- data.frame(accvec,sizes, cluster=names(accvec))
df <- df[order(df$sizes),]


ggplot(df, aes(sizes,accvec, label=cluster))+
  geom_point()+
  ggrepel::geom_text_repel(size = 1.5)+
  xlab('cluster size') + ylab('accuracy')


