library(tsne)
library(randomcoloR)
source("genovese_scripts/modes.functions.r")

normalize<-function(m){
  # Normalizes a matrix my column
  (m - mean(m, na.rm = TRUE))/sd(m, na.rm = TRUE)
}

normByMouse<-function(data.sc, featID){
  dataNorm = data.sc
  for (M in unique(data.sc$id)){
    mdat = data.sc[data.sc$id == M, featID]
    dataNorm[data.sc$id == M,1:ncol(mdat)] <- as.data.frame(lapply(mdat, normalize))
  }
  return(dataNorm)
}


KDECluster<-function(X,NIter,BW){
  # Computes the optimal number of density modes in an abstract feature space.
  #  This function clusters an arbitray dataspace using the mean shift algorithm on a kernel density manifold.
  #  You may pass a spectrum of bandwidths to assess the dataset for cluster number and return...
  #  cluster IDs for each datapoint at the critical bandwidth (Bandwidth with most sig modes).
  #  OR you may pass a single bandwidth to mean shift the dataset at that bandwidth. The function will return
  #  Cluster IDs and the significance of modes identified but eb unable to comment on the critical width.
  #  
  #
  # Args:
  #   X: Must be MATRIX of observations x variables. Variables may be in any feature space.
  #   NIter: Number of iterations of the clustering algorithm. DEFAULT = 1. TO DO
  #   BW: Spectrum of bandwidths to test KDE modes over.
  #
  # Returns:
  #   stats: A list of KDE clustering features including:
  #     optimalBW: Bandwidth at which the greatest number of significant modes are...
  #                present of those tested.
  #     ClusterID: An array of integers representing cluster assignment at ...
  #                 the optimal bandwidth of in BW.
  #     NModes: An array of integers with the number of modes in datX at each bandwidth.
  #     NSigModes: An array of integers with number of modes found to be significant at each bandwidth.
  #     SigModeLoc: The location in featuer space of the significant modes found in datX.
  #     NonSigModeLoc: The location in featuer space of the non-significant modes found in datX.
  #     AggMode: The location in featuer space of the non-significant modes found in datX.
  # 
  # To Do: OSS Alter datastructures to facilitate NIter > 1
  
  
  for(I in 1:NIter){
    samp = sample(x=1:nrow(X), size=nrow(X), replace=FALSE)
    datX = X[head(samp,nrow(X)/2),1:ncol(X)]
    datY = X[head(samp,nrow(X)/2),1:ncol(X)]
    
    SigModes = matrix(,nrow=100,ncol=ncol(X))
    NSigModes = matrix(,nrow=length(BW),ncol=1)
    NTotalModes = matrix(,nrow=length(BW),ncol=1)
    ClusterID = matrix(,nrow=nrow(X),ncol=1)
    
    print('Computing bandwidth sweep, this may take some time depending on the number of modes...')
    XMBT = multibandtest.fun(datX, datY, BW, nboot=1000, alpha=0.05, digits=2)
  
    # Extract number of significant modes per bandwidth..
    for(h in 1:length(BW)){
  	  NSigModes[h] = sum(XMBT[[h]][["CIlead"]][,2] < 0)
  	  NTotalModes[h] = nrow(XMBT[[h]][["CIlead"]])	
    }
    
    print(paste('Maximum no. of significant modes found = ',max(NSigModes),sep=''))
    bwidx = which(!is.na(match(NSigModes, max(NSigModes)))) # Identify Critical bandwidth as that with most significant modes
    bwidx = bwidx[1] # 1 For smallest length(bwidx) for greatest
    sidx = XMBT[[bwidx]][["CIlead"]][,2] < 0 # Find indices of significant modes
    nsidx = XMBT[[bwidx]][["CIlead"]][,2] > 0 # Find indices of non-significant modes
    SigModeLoc = matrix(XMBT[[bwidx]][['modes']][sidx,],ncol=ncol(X)) # Find significant mode locations
    NonSigModeLoc = matrix(XMBT[[bwidx]][['modes']][nsidx,],ncol=ncol(X)) # Find non-significant mode locations
    
    print(paste('Assigning cluster identities at critical bandwidth, may take time...',sep=''))
    critical = round(matrix(t(apply(X=X, MARGIN=1, FUN=msiter.fun, dat=X, bw=BW[bwidx])), ncol=ncol(X)),digits=2) # Finds the cluster origin for every data point at the critical bandwidth
    AggModes = modecl.fun(modemat=critical, digits=2)
    print(paste('Completed...',sep=''))
    
    # Assign Cluster ID based on modes found in entire dataset
    for(pt in 1:nrow(critical)){
  	  for(M in 1:nrow(AggModes)){
  		  if(sum(AggModes[M,] == critical[pt,])==ncol(X)){
  			  ClusterID[pt] = M;
  		  }
  	  }
    }
  }
  stats <- list(optimalBW = BW[bwidx], ClusterID = ClusterID, NTotalModes = NTotalModes, NSigModes = NSigModes, SigModeLoc = SigModeLoc, NonSigModeLoc = NonSigModeLoc, AggModes = AggModes)
  return(stats)
}


KDESummaryFig<-function(data.sc,PCstats,TSNEstats,dr){
  # Summary FIgure Plotting
  TSNEplotCol = distinctColorPalette(k = nrow(TSNEstats$AggModes), altCol = TRUE, runTsne = FALSE)
  PCplotCol = distinctColorPalette(k = nrow(PCstats$AggModes), altCol = TRUE, runTsne = FALSE)
  expCol = distinctColorPalette(k = length(unique(data.sc.norm$expr)), altCol = FALSE, runTsne = FALSE)
  mouseCol = distinctColorPalette(k = length(unique(data.sc.norm$id)), altCol = TRUE, runTsne = FALSE)
  
  
  psize = 1.2
  msize = 2.5
  
  pdf(file=paste('figures/','SCvsPCsummaryFig.pdf',sep=''),
      width=25, height=10)
  par(mfrow=c(2,5))
  # Plot 1
  plot(c(),main='t-SNE Map',xlab='t-SNE 1',ylab='t-SNE 2',xlim=range(dr[,1])*1.2,ylim=range(dr[,2])*1.2)
  points(dr[,1],dr[,2],pch=20,cex=psize,col='black')
  
  # Plot 2
  plot(c(),main='t-SNE by t-SNE Cluster ID',xlab='t-SNE 1',ylab='t-SNE 2',xlim=range(dr[,1])*1.2,ylim=range(dr[,2])*1.2)
  for(C in unique(TSNEstats$ClusterID)){
    points(dr[TSNEstats$ClusterID==C,1],dr[TSNEstats$ClusterID==C,2],pch=20,cex=psize,col=TSNEplotCol[C]) # Plot t-sne points
  }
  points(TSNEstats$SigModeLoc[,1],TSNEstats$SigModeLoc[,2],pch=20,cex=msize,col='dodgerblue4') # Plot non-sig modes
  points(TSNEstats$NonSigModeLoc[,1],TSNEstats$NonSigModeLoc[,2],pch=20,cex=msize,col='firebrick2') # Plot Sig modes
  points(TSNEstats$AggModes[,1],TSNEstats$AggModes[,2],pch=20,cex=msize,col='black') # Plot full dataset modes
  text(TSNEstats$AggModes[,1]+0.05*diff(range(dr[,1])),TSNEstats$AggModes[,2]+0.05*diff(range(dr[,1])), seq(1,nrow(TSNEstats$AggModes),length=nrow(TSNEstats$AggModes)),cex=2) # Plot labels
  
  # Plot 3
  plot(c(),main='t-SNE by PC Cluster ID',xlab='t-SNE 1',ylab='t-SNE 2',xlim=range(dr[,1])*1.2,ylim=range(dr[,2])*1.2)
  for(C in unique(PCstats$ClusterID)){
    points(dr[PCstats$ClusterID==C,1],dr[PCstats$ClusterID==C,2],pch=20,cex=psize,col=PCplotCol[C]) # Plot t-sne points
  }
  
  # Plot 4
  plot(c(),main='t-SNE by Experimenter',xlab='t-SNE 1',ylab='t-SNE 2',xlim=range(dr[,1])*1.2,ylim=range(dr[,2])*1.2)
  i=1
  for(E in unique(data.sc.norm$expr)){
    points(dr[data.sc.norm$expr==E,1],dr[data.sc.norm$expr==E,2],pch=20,cex=psize,col=expCol[i]) # Plot t-sne points
    i=i+1
  }
  legend("topright",legend=c("HP", "DG"),col=expCol,pch=19,cex=1.4)
  
  # Plot 5
  plot(c(),main='t-SNE by Celltype',xlab='t-SNE 1',ylab='t-SNE 2',xlim=range(dr[,1])*1.2,ylim=range(dr[,2])*1.2)
  i=1
  for(M in unique(data.sc.norm$id)){
    points(dr[data.sc.norm$id==M,1],dr[data.sc.norm$id==M,2],pch=20,cex=psize,col=mouseCol[i]) # Plot t-sne points
    i=i+1
  }
  legend("topright",legend=c("Calbindin", "SC"),col=mouseCol,pch=19,cex=1.4) # Relevant only for cell type comparisons
  
  # Plot 6
  plot(c(),main='PC Map',xlab='PC 1',ylab='PC 2',xlim=range(data.sc.pca[,1])*1.2,ylim=range(data.sc.pca[,2])*1.2)
  points(data.sc.pca[,1],data.sc.pca[,2],pch=20,cex=psize,col='black')
  
  # Plot 7
  plot(c(),main='PC by t-SNE Cluster ID',xlab='PC 1',ylab='PC 2',xlim=range(data.sc.pca[,1])*1.2,ylim=range(data.sc.pca[,2])*1.2)
  for(C in unique(TSNEstats$ClusterID)){
    points(data.sc.pca[TSNEstats$ClusterID==C,1],data.sc.pca[TSNEstats$ClusterID==C,2],pch=20,cex=psize,col=TSNEplotCol[C]) # Plot PC points
  }
  
  # Plot 8
  plot(c(),main='PC by PC Cluster ID',xlab='PC 1',ylab='PC 2',xlim=range(data.sc.pca[,1])*1.2,ylim=range(data.sc.pca[,2])*1.2)
  for(C in unique(PCstats$ClusterID)){
    points(data.sc.pca[PCstats$ClusterID==C,1],data.sc.pca[PCstats$ClusterID==C,2],pch=20,cex=psize,col=PCplotCol[C]) # Plot PC points
  }
  points(PCstats$SigModeLoc[,1],PCstats$SigModeLoc[,2],pch=20,cex=msize,col='dodgerblue4') # Plot non-sig modes
  points(PCstats$NonSigModeLoc[,1],PCstats$NonSigModeLoc[,2],pch=20,cex=msize,col='firebrick2') # Plot Sig modes
  points(PCstats$AggModes[,1],PCstats$AggModes[,2],pch=20,cex=msize,col='black') # Plot full dataset modes
  text(PCstats$AggModes[,1]+0.02*diff(range(data.sc.pca[,1])),PCstats$AggModes[,2]+0.02*diff(range(data.sc.pca[,1])), seq(1,nrow(PCstats$AggModes),length=nrow(PCstats$AggModes)),cex=2) # Plot labels
  
  # Plot 9
  plot(c(),main='PC by Experimenter',xlab='PC 1',ylab='PC 2',xlim=range(data.sc.pca[,1])*1.2,ylim=range(data.sc.pca[,2])*1.2)
  i=1
  for(E in unique(data.sc.norm$expr)){
    points(data.sc.pca[data.sc.norm$expr==E,1],data.sc.pca[data.sc.norm$expr==E,2],pch=20,cex=psize,col=expCol[i]) # Plot t-sne points
    i=i+1
  }
  legend("topright",legend=c("HP", "DG"),col=expCol,pch=19,cex=1.4)
  
  # Plot 10
  plot(c(),main='PC by Cell Type',xlab='PC 1',ylab='PC 2',xlim=range(data.sc.pca[,1])*1.2,ylim=range(data.sc.pca[,2])*1.2)
  i=1
  for(M in unique(data.sc.norm$id)){
    points(data.sc.pca[data.sc.norm$id==M,1],data.sc.pca[data.sc.norm$id==M,2],pch=20,cex=psize,col=mouseCol[i]) # Plot t-sne points
    i=i+1
  }
  legend("topright",legend=c("Calbindin", "SC"),col=mouseCol,pch=19,cex=1.4) # Relevant only for cell type comparisons
  dev.off()
}

