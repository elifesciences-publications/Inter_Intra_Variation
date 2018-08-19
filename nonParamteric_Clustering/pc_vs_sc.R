source("KDEfunc.R")
library(tidyverse)
library(tsne)

# Import PC and SC cell feature datafiles
fname.sc <- "raw_data/SC_OSSdatatable.txt" # Import NOLANLab .txt datafile
fname.pc <- "raw_data/PC_OSSdatatable.txt"

data.scimport <- read_tsv(fname.sc)
data.sc <- data.scimport %>% drop_na(dvloc)

data.pcimport <- read_tsv(fname.pc)
data.pc <- data.pcimport %>% drop_na(dvloc)

# Convert DV location to mm
data.sc <- mutate(data.sc, dvlocmm = dvloc/1000)
data.pc <- mutate(data.pc, dvlocmm = dvloc/1000)

# Assign cell type labels
data.sc$celltype <- rep(0,nrow(data.sc))
data.pc$celltype <- rep(1,nrow(data.pc))

# Downsample SCs and form combined dataset
samp = sample(1:nrow(data.pc), nrow(data.pc))
data.comb <- rbind(data.sc[samp,1:11], size = data.pc[1:11])

# Normalise data and add label vars
data.comb.norm <- as.data.frame(lapply(data.comb[1:11], normalize)) # Normalize by all
data.comb.norm$dvlocmm <- c(data.sc$dvlocmm[samp], data.pc$dvlocmm)
data.comb.norm$id <- c(data.sc$mid[samp], data.pc$mid)
data.comb.norm$expr <- c(data.sc$expr[samp], data.pc$expr)
data.comb.norm$id <- c(data.sc$celltype[samp], data.pc$celltype)


# Compute principal components
pca <- prcomp(data.comb.norm[1:11])
data.comb.pca <- as.data.frame(pca[["x"]])

# Assign dataspace to be used in t-SNE and clustering analysis
dat = data.comb.pca[1:3]
dat = as.matrix(dat,nrow(dat),ncol(dat))

# Remap into 2D t-SNE Space
print('Computing t-SNE map')
dr = tsne(dat, initial_config = NULL, k = 2, initial_dims = ncol(dat), perplexity = 20, max_iter = 1000, min_cost = 0, epoch_callback=NULL, whiten = TRUE, epoch=100)

# KDE Cluster in 2D t-SNE Space
print('Clustering in 2D t-SNE space')
TSNEstats = KDECluster(X=dr,NIter=1,BW=seq(80, 90, length=10))
# bw = 4.7 for sc 3 cluster form

# KDE Cluster in 3D principal component space
print('Clustering in 3D PC space')
PCstats = KDECluster(X=dat,NIter=1,BW=seq(0.9, 1.5, length=10))
# bw = 0.9 for sc 3 cluster form

# Assign t-SNE identified cluster IDs to dataframe for export
data.comb$OSSClusterID <- TSNEstats$ClusterID
write.table(data.comb,file="raw_data/COMBdatatablewClusters.txt",sep="\t",row.names=FALSE) # Export as .txt file

# Reassign combined data to 'sc' label for use in generic summary fitting function...
data.sc = data.comb
data.sc.norm = data.comb.norm
data.sc.pca = data.comb.pca

# Construct summary fig with remappings and cluster IDs
KDESummaryFig(data.sc,PCstats,TSNEstats,dr)

