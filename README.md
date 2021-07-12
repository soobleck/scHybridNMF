# scHybridNMF
single-cell hybrid NMF package used to bicluster gene expression and spatial location data of single cells.

Files include:
- **scHybridNMF**: a combined NMF scheme that takes:
	-- a 2-by-n location matrix, 
	-- an m-by-n gene expression matrix of single cell data, 
	-- the number of clusters k,
	-- (OPTIONAL) a location-gene expression weighting term alpha,
	-- (OPTIONAL) the tolerance level for the normalized KKT residual stopping criterion,
	-- (OPTIONAL) a maximum number of iterations.
- All other files are for running NMF.