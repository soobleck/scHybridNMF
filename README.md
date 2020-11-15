# scHybridNMF
single-cell hybrid NMF package used to bicluster gene expression and spatial location data of single cells.

Files include:
- **scHybridNMF**: a combined NMF scheme that takes a pair of 2-by-n location matrix and m-by-n gene expression matrix of single cell data.
- **vis_clusters**: a script to save side-by-side dot plots of cell locations and t-SNE plots of gene expression data. The dots are color-coded based on cluster membership.
- All other files are for running NMF.