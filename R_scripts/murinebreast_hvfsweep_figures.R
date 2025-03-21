library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(ComplexHeatmap)
library(grid)


set.seed(2022)




#read in murine breast sweep result table from Rohan



### 4B
dat_4b <- read.csv('Files for Alex/figure_4b_violin_data.csv')
dat <- dat_4b

dat <- dat[order(dat$leiden),]


#prep for plot
dat <- reshape2::melt(dat, id.vars = 'leiden')

#remove NAs
dat <- dat[complete.cases(dat),]

#prop to percent
dat$value <- dat$value * 100

#remove prefix "X"
dat$features <- plyr::mapvalues(dat$variable,
                                from = levels(dat$variable),
                                to = as.character(gsub(x=levels(dat$variable), 'X', '')))



g <- ggplot(dat, aes(features, value))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'steelblue', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  theme_light()+
  ylab('% accuracy') + xlab('Binary LSI Clusters')+
  theme(axis.text = element_text(size = 8),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )



pdf('Files for Alex/4b_new_vln.pdf', height =3, width = 4)

g
# g +scale_y_continuous(limits = c(20,100))

dev.off()




## 4f - rna sweep 
dat_4f <- read.csv('Files for Alex/breast_cluster_accuracy_violin_data_rna_varying_mofa_figure4F.csv')
dat <- dat_4f

#prep for plot
dat <- reshape2::melt(dat, id.vars = 'modified_hvf_lsi_cluster')

dat$features <- plyr::mapvalues(dat$variable,
                                from = levels(dat$variable),
                                to = c("0", "800", '1500', '2300', '3100', '2800', '4600', '5400', '6100', '6900', '7700'))

g <- ggplot(dat, aes(features, value))+
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



pdf('Files for Alex/4f_new_RNA.pdf', height =3, width = 4)

g
g +scale_y_continuous(limits = c(20,100))

dev.off()





## 4g - atac sweep 
dat_4g <- read.csv('Files for Alex/breast_cluster_accuracy_violin_data_atac_varying_mofa_figure4G.csv')
dat <- dat_4g

#prep for plot
dat <- reshape2::melt(dat, id.vars = 'modified_hvf_lsi_cluster')

dat$features <- plyr::mapvalues(dat$variable,
                                from = levels(dat$variable),
                                to = c("0", '1200', '2500', '3700', '4900', '6200', '7400', '8600', '9900', '11100', '12300'))

g <- ggplot(dat, aes(features, value))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'firebrick', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  theme_light()+
  ylab('Cluster Dominant Celltype Accuracy') + xlab('Number of included ATAC features sorted by variance')+
  theme(axis.text = element_text(size = 5),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )+
  ggtitle('Varying ATAC features, effect on Cluster Accuracy')



pdf('Files for Alex/4g_new_ATAC.pdf', height =3, width = 4)

g
g +scale_y_continuous(limits = c(20,100))

dev.off()






## s8 atac
dat_s8 <- read.csv('Files for Alex/fig_S8_breast_cluster_accuracy_violin_data_atac_varying_figure_vs_integ.csv')
dat <- dat_s8

# dat <- dat[order(dat$modified_hvf_lsi_cluster),]


#prep for plot
dat <- reshape2::melt(dat, id.vars = 'modified_hvf_lsi_cluster')

#remove NAs
# dat <- dat[complete.cases(dat),]

#prop to percent
# dat$value <- dat$value * 100

#remove prefix "X"
dat$features <- plyr::mapvalues(dat$variable,
                                from = levels(dat$variable),
                                to = as.character(gsub(x=levels(dat$variable), 'X', '')))



g <- ggplot(dat, aes(features, value))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'firebrick', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  theme_light()+
  ylab('% accuracy of clusters') + xlab('All RNA HVF + number of ATAC peak HVFs used')+
  ggtitle('Cluster Accuracy by number of ATAC Peaks used')+
  theme(axis.text = element_text(size = 7),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )



pdf('Files for Alex/s8_atac.pdf', height =3, width = 4)

g
# g +scale_y_continuous(limits = c(20,100))

dev.off()





## s8 rna
dat_s8 <- read.csv('Files for Alex/fig_S8_breast_cluster_accuracy_violin_data_rna_varying_figure_vs_integ.csv')
dat <- dat_s8

# dat <- dat[order(dat$modified_hvf_lsi_cluster),]


#prep for plot
dat <- reshape2::melt(dat, id.vars = 'modified_hvf_lsi_cluster')

#remove NAs
# dat <- dat[complete.cases(dat),]

#prop to percent
# dat$value <- dat$value * 100

#remove prefix "X"
dat$features <- plyr::mapvalues(dat$variable,
                                from = levels(dat$variable),
                                to = as.character(gsub(x=levels(dat$variable), 'X', '')))



g <- ggplot(dat, aes(features, value))+
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.7)+
  geom_violin(fill = 'steelblue', alpha = 0.3, linewidth = 0.5)+
  geom_boxplot(alpha = 0, width = 0.2)+
  ggtitle('Cluster Accuracy by number of RNA Genes used')+
  theme_light()+
  ylab('% accuracy of clusters') + xlab('All ATAC HVF + number of RNA peak HVFs used')+
  theme(axis.text = element_text(size = 7),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(face = "bold", size = 10)
  )



pdf('Files for Alex/s8_rna.pdf', height =3, width = 4)

g
# g +scale_y_continuous(limits = c(20,100))

dev.off()
