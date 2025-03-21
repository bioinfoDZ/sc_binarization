library(tidyverse)
library(Seurat)
library(msigdbr)
library(ComplexHeatmap)
library(RColorBrewer)
library(foreach)
library(parallel)
library(Matrix)




### raw data downloaded from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128423

### scbroad, download metadata and some normalized data
# https://singlecell.broadinstitute.org/single_cell/study/SCP361/mouse-bone-marrow-stroma-in-homeostasis#study-visualize


#read in tsne and cluster info from paper
tsne <- data.table::fread('data/bmstroma/broadscportal/SCP361/cluster/stroma.tsne.txt', header =T)
tsne <- tsne[-1,]
colnames(tsne)[2:3] <- c('tsne_1', 'tsne_2')

clust <- data.table::fread('data/bmstroma/broadscportal/SCP361/metadata/stroma.tsne.meta.txt', header = T)
clust <- clust[-1,]


#check if order matching
table(clust$NAME == tsne$NAME)
md <- cbind(clust, tsne[,-1])

rm(clust,tsne)

#make dataframe, just in case...
md <- as.data.frame(md)

#Name -> barcode
colnames(md)[1] <- 'barcode'

#get sample name info from the md
md$Sample <- str_split_fixed(md$barcode, '_', 2)[,1]

#change order, not really important but just want sample ahead
md <- md[,c(1,6,2,3,4,5)]

#change object types

#sample and clusters: make factors, levels by number
md[,2] <- factor(md$Sample)
md[,3] <- factor(md$Cluster, levels = names(sort(table(md$Cluster), decreasing = T) ))
md[,4] <- factor(md$Subcluster, levels = names(sort(table(md$Subcluster), decreasing = T) ))

#tsne embeddings, numeric
md[,5] <- as.numeric(md[,5])
md[,6] <- as.numeric(md[,6])



#add in cell types from paper
celltype <- data.frame(cluster = levels(md$Cluster),
                       celltype = c('EC-Sinusoidal', 'Lepr-MSC', 'Chondro-hyper', 'Fibro-4', 'Chondro-progen', 'Fibro-5',
                                    'EC-arteriolar', 'OLC-1', 'OLC-2', 'Fibro-1', 'Chondro-prol/rest',
                                    'EC-arterial', 'Pericytes', 'Chondro', 'Fibro-2', 'Fibro-3', 'Chondro-prehyper-2'
                                    
                       )
                       
)

#add in anatomic classification from paper...
celltype$anat_class <- 'Bone'
celltype[celltype$cluster %in% c(1,7,12,0,6,11),"anat_class"] <- 'BM'

#add celltype to md
md$Celltype <- md$Cluster
md$Celltype <- plyr::mapvalues(md$Celltype, levels(md$Cluster), celltype$celltype)

#add anatomic classification to md
md$anat_class <- md$Cluster
md$anat_class <- plyr::mapvalues(md$anat_class, levels(md$Cluster), celltype$anat_class)


#plot tsne
ggplot(md, aes(tsne_1,tsne_2, col = Sample))+
  geom_point(size = 0.1)


#make repel labels...
repel <- aggregate(tsne_1 ~ Cluster, md, median)
repel$tsne_2 <- aggregate(tsne_2 ~ Cluster, md, median)[,2]
ggplot(md, aes(tsne_1,tsne_2, col = Cluster))+
  geom_point(size=0.1)+
  ggrepel::geom_text_repel(data = repel, inherit.aes = F, aes(tsne_1, tsne_2, label = Cluster))

ggplot(md, aes(tsne_1,tsne_2, col = Subcluster))+
  geom_point(size=0.1)


ggplot(md, aes(tsne_1,tsne_2, col = anat_class))+
  geom_point(size=0.1)

#make repel labels...
repel <- aggregate(tsne_1 ~ Celltype, md, median)
repel$tsne_2 <- aggregate(tsne_2 ~ Celltype, md, median)[,2]
ggplot(md, aes(tsne_1,tsne_2, col = Celltype))+
  geom_point(size=0.1)+
  ggrepel::geom_text_repel(data = repel, inherit.aes = F, aes(tsne_1, tsne_2, label = Celltype))




#######  read in the data ######
# GEO is missing two samples??? std7, std8???
# scportal has "normalized" data
# it is "TP10K" normalized
# Gene expression was represented as the fraction of its UMI count 
# with respect to total UMI in the cell and then multiplied by 10,000. 
# We denoted it by TP10K – transcripts per 10K transcripts.


# this is similar to "RC" normalization method.
# log1p on top of this should recapitulate normal seurat log normalization from raw counts.

#read in normalized
dat <- data.table::fread('data/bmstroma/broadscportal/SCP361/expression/stroma.TP4K.txt')

#use genes as rownames, remove gene col
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]



#read as seurat object
# Counts are not really counts here, they are log1p relative counts...
sobj <- CreateSeuratObject(dat, project = 'bmstroma')
sobj@assays$RNA@data <- log1p(sobj@assays$RNA@counts)
rm(dat)

#add in metadata
sobj@meta.data <- cbind(sobj@meta.data, md[,-1])

#set default idents to celltypes
sobj <- SetIdent(sobj, value = sobj$Celltype)


### add in tsne embeddings
# get embeddings, set rownames = barcodes, set as matrix, turn to seurat dimreduc object, then add to sobj
tsne <- md[,grepl('tsne', colnames(md))]
rownames(tsne) <- colnames(sobj)
tsne <- as(tsne, 'matrix')
tsne <- Seurat::CreateDimReducObject(embeddings = tsne, global=T)
sobj[['tsne']] <- tsne



#save it
saveRDS(sobj, 'data/bmstroma/parsed_sobj.rds')





