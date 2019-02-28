source("nonParametric_Clustering/KDEfunc.R")
library(tidyverse)
library(tsne)

# Import cell feature datafile
fname.sc <- "nonParametric_Clustering/raw_data/SC_HPdatatable.txt" # Import NOLANLab .txt datafile
data.import <- read_tsv(fname.sc)
data.sc <- data.import %>% drop_na(dvloc)

# Normalise as a whole and by mouse
data.sc <- mutate(data.sc, dvlocmm = dvloc/1000)
data.sc.norm <- as.data.frame(lapply(data.sc[1:11], normalize)) # Normalize by all
data.sc.Mnorm <- normByMouse(data.sc, 1:11)
data.sc.norm$dvlocmm <- data.sc$dvlocmm
data.sc.norm$id <- data.sc$id
data.sc.norm$housing <- data.sc$housing
data.sc.norm$expr <- data.sc$expr

# Compute principal components
pca <- prcomp(data.sc.norm[1:11])
data.sc.pca <- as.data.frame(pca[["x"]])

# Assign dataspace to be used in t-SNE and clustering analysis
dat = data.sc.pca[1:3]
dat = as.matrix(dat,nrow(dat),ncol(dat))

# Remap into 2D t-SNE Space
print('Computing t-SNE map')
dr = tsne(dat, initial_config = NULL, k = 2, initial_dims = ncol(dat), perplexity = 80, max_iter = 1000, min_cost = 0, epoch_callback=NULL, whiten = TRUE, epoch=100)

# KDE Cluster in t-SNE Space
print('Clustering in 2D t-SNE space')
TSNEstats = KDECluster(X=dr,NIter=1,BW=seq(4.7, 4.7, length=1))

# KDE Cluster in Original principal component space
print('Clustering in 3D PC space')
PCstats = KDECluster(X=dat,NIter=1,BW=seq(0.9, 0.9, length=1))

# Assign t-SNE identified cluster IDs to dataframe for export
data.sc$OSSClusterID <- TSNEstats$ClusterID
write.table(data.sc,file="raw_data/datatablewClusters.txt",sep="\t",row.names=FALSE) # Export as .txt file

# Construct summary fig with remappings and cluster IDs
KDEsum <- KDESummaryFig(data.sc,PCstats,TSNEstats,dr)
