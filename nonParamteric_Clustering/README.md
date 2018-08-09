## Synopsis
This is a repository for the Kernel Density Estimate Clustering undertaken on a large patch clamp dataset from stellate cells in mouse medial entorhinal cortex. The main script is KDEClustering.Rmd and explains the step-by-step procedure and plots relevant figures.

## Directory Structure

### README.md

### KDEClustering.Rmd
	Central script narratively describing the method and relevant code.
### cell_att_scripts
	Matlab scripts required for extracting features from raw electrophys traces in .h5 	format. Original extraction scripts were written by Hugh Pastoll.
### raw_data
	.txt files with cell features and descriptives in. datatable.txt is a smaller more refined selection of stellate cells. OSSFeatDat is features of all cells originally patched. those marked wClusters have an extra variable named OSScluster which are the cluster identities from KDE clustreing undertaken on the larger OSSFeatDat cells.
### figures
	.pdf figures of custer results on both data sets.
### genovese_scripts
	Scripts written by Genovese et al. in their paper "Non-	parametric inference for 		density modes":
	Genovese, C.R., Perone-Pacifico, M., Verdinelli, I., Wasserman, L., 2015. Non-parametric inference for density modes. J. R. Stat. Soc. Ser. B Stat. Methodol. https://doi.org/10.1111/rssb.12111
	Availble here: [link] https://sites.google.com/a/uniroma1.it/marcoperonepacifico/R-code.