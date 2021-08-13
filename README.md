# scHybridNMF
single-cell hybrid NMF package used to bicluster gene expression and spatial location data of single cells.

Files include:
- **scHybridNMF**: a combined NMF scheme that takes:
	1. a 2-by-n location matrix, 
	2. an m-by-n gene expression matrix of single cell data, 
	3. the number of clusters k,
	4. (OPTIONAL) a location-gene expression weighting term alpha,
	5. (OPTIONAL) the tolerance level for the normalized KKT residual stopping criterion,
	6. (OPTIONAL) a maximum number of iterations.
- All other files are for running NMF.
