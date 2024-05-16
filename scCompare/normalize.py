import numpy as np
import pandas as pd

import tqdm

import pyarrow

pyarrow.set_io_thread_count(20)#FIXME this is not working as I want it to (+remove hard-coded)

def load_gene_cell_expr(expr_mat, remove_zeros=True):
	"""
	Load gene-cell expression matrix (cells are columns, rows genes).
	"""
	print('--- Loading gene expression per cell ---')
	expr = pd.read_csv(expr_mat, sep="\t", engine='pyarrow', header=0, index_col=0)

	#Filter out non-expressed genes
	# expr = expr.loc[~(expr==0).all(axis=1)]

	#Filter out lowly-expressed genes (at least 10 UMI)
	expr = expr.loc[(expr.sum(axis=1) >= 10)]
	return expr


def load_cell_clust(cells_to_clusters, filter_out_start=None):
	"""
	Load cell to clusters
	"""
	clust = pd.read_csv(cells_to_clusters, sep="\t", header=0, index_col=0)
	if filter_out_start:
		clust = clust[~clust['cluster_name'].astype(str).str.startswith(filter_out_start)]
	clust = clust['cluster_name'].to_dict()
	val = sorted(set(clust.values()))
	return clust, val


def geometric_mean(vector):
	return (np.exp(np.mean(np.log(1 + vector)))) - 1


def normalize_geom_mean(expr, clust, val):

	print('--- Computing normalized gene expression per cluster ---')

	# For each cluster, compute gene expression as geometric mean over cells (do parallel?)
	data = []
	tot = 0
	for v in tqdm.tqdm(val, colour='#FF79C6'):
		cells = [i for i in clust if clust[i]==v]
		lg = len(cells)
		tot += lg
		# print(f'Cluster {v}: {lg} cells')
		expr_cells = expr[cells]

		mat = np.array(expr_cells)
		mc_meansize = np.mean(mat.sum(axis=0))
		ideal_cell_size = min(1000, np.median(mc_meansize))
		# print(ideal_cell_size)
		gene_expr_in_clust = np.zeros(len(mat[:,0]))
		for i in range(len(mat[:,0])): #if we keep the code, remove the loop and do matrix operation (operation on each row)
			v = mat[i,:]
			regmean = (geometric_mean(v) / mc_meansize) * ideal_cell_size
			gene_expr_in_clust[i] = regmean
		data.append(gene_expr_in_clust)


	# Normalize gene expression by dividing by its median of expression across clusters
	data = np.array(data).T
	medians = np.median(data, axis=1)
	for i in range(len(medians)):
		data[i,:] = (data[i,:] + 0.05) / (medians[i] + 0.05) #if we keep the code, remove the loop and do matrix operation
	
	print(f'Total: {tot} cells')
	return data


def save_outputs(data, genes, clusters, output):
	df = pd.DataFrame(data=data, index=genes, columns=clusters)
	df.to_csv(output, sep='\t')


def normalize_main(expr_mat, cells_to_clusters, output, filter_out_start=None):
	matrix = load_gene_cell_expr(expr_mat)
	clusters, cluster_names = load_cell_clust(cells_to_clusters, filter_out_start=filter_out_start)
	norm_matrix = normalize_geom_mean(matrix, clusters, cluster_names)
	save_outputs(norm_matrix, matrix.index, cluster_names, output)


# Metacell normalization
# sca_cell_type_fp <- function( input_table, mc_object, mat_object, nbins=10L) {
	
# 	# load input table
# 	if (any(class(input_table) %in% "character")) {
		
# 		cell_type_table = read.table(input_table, header = TRUE, sep="\t", comment.char="")
# 		colnames(cell_type_table) = c("metacell", "cell_type", "color")
# 		cell_types_ordered = unique(cell_type_table$cell_type)
# 		rownames(cell_type_table) = cell_type_table$metacell
		
# 	} else if (any(class(input_table) %in% "data.frame")) {
		
# 		cell_type_table = input_table
# 		class(cell_type_table) <- "data.frame"
# 		colnames(cell_type_table) = c("metacell", "cell_type", "color")
# 		cell_types_ordered = unique(cell_type_table$cell_type)
# 		rownames(cell_type_table) = cell_type_table$metacell
		
# 	}
	
# 	sc_ct_label = as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
# 	names(sc_ct_label) = names(mc_object@mc)
	
# 	cells_cols = cell_type_table[as.character(mc_object@mc),"color"]
# 	cells_cols = as.character(cells_cols)
# 	names(cells_cols) = names(mc_object@mc)
	
# 	# filter low expression genes not included in the mc_fp
# 	umis = mat_object@mat
# 	umis = umis[rownames(mc_object@mc_fp),names(mc_object@mc)]
	
# 	#ct_counts=t(apply(umis,1,function(x) tapply(x, sc_ct_label, sum)))
# 	#ct_size=colSums(ct_counts)
# 	#ct_umifrac=t(apply(ct_counts,1,function(x) x*1000/ct_size))
# 	#ct_umifrac_n=(0.1 + ct_umifrac)/apply(0.1 + ct_umifrac, 1, median)
	
# 	ct_geomean = tryCatch(
# 	  t(apply(umis, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1))),
# 	  error = function(e) {
# 	    warning(e)
# 	    message("Calculating geom mean for genes in ",nbins," bins")
# 	    umis_list <- vector("list",nbins)
# 	    sl <- split(1:nrow(umis),cut(1:nrow(umis),nbins))
# 	    for (i in seq_along(sl)) {
# 	      message(i," / ", nbins)
# 	      umisub <- umis[sl[[i]],]
# 	      ct_geomean_sub <- t(apply(umisub, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1)))
# 	      umis_list[[i]] <- ct_geomean_sub 
# 	    }
# 	    do.call(rbind, umis_list)
# 	  }
# 	)
# 	ct_geomean = ct_geomean[ , cell_types_ordered ]
# 	ct_meansize = tapply(Matrix::colSums(umis), sc_ct_label, mean)
# 	ideal_cell_size = pmin(1000,median(ct_meansize))
# 	g_fp = t(ideal_cell_size * t(ct_geomean) / as.vector(ct_meansize))
# 	fp_reg = 0.05
# 	g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
	
# 	# for compatibility with other functions, return as a MC-like object
# 	ct_table=mc_object
# 	ct_table@mc_fp=g_fp_n
# 	ct_table@mc=sc_ct_label
# 	ct_table@colors=cells_cols
# 	ct_table@cell_names=names(sc_ct_label)
# 	return(ct_table)
	
# }



# #' Compute metacell gene footprint
# #'
# #' The footprint is defined as the size-normalized geometric mena of the number of umis per metacells, dvided by the median over all metacells.
# #'
# #' @param mc a metacell object
# #' @param us umi matrix
# #' @param norm_by_mc_size normalize by mean total umis and then multiply by median mean mc size (or at least 1000). This means umis per 1000 molecules.
# #' @param min_total_umi consider genes with at least min_total_umi total umis
# #'
# #' @export
# mc_compute_fp = function(mc, us, norm_by_mc_size=T, min_total_umi=10)
# {
# 	f_g_cov = rowSums(us) > min_total_umi

# 	if(0) {
# 		mc_cores = get_param("mc_cores")
# 		doMC::registerDoMC(mc_cores)
# 		all_gs = rownames(us[f_g_cov,])
# 		n_g = length(all_gs)
# 		g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
# 		fnc = function(gs) {
# 						.row_stats_by_factor(us[gs,],
# 										mc@mc,
# 										function(y) {exp(rowMeans(log(1+y)))-1}) }

# 		clust_geomean = do.call(rbind, mclapply(g_splts, fnc, mc.cores = mc_cores))
# 	}
# 	clust_geomean = t(tgs_matrix_tapply(us[f_g_cov,], mc@mc,
# 									  function(y) {exp(mean(log(1+y)))-1}))
# 	rownames(clust_geomean) = rownames(us)[f_g_cov]

# #	clust_geomean = .row_stats_by_factor(us[f_g_cov,],
# #									mc@mc,
# #									function(y) {exp(rowMeans(log(1+y)))-1})

# 	if (norm_by_mc_size) {
# 		mc_meansize = tapply(colSums(us), mc@mc, mean)
# 		ideal_cell_size = pmin(1000, median(mc_meansize))
# 		g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
# 	}
# 	else {
# 		g_fp = clust_geomean
# 	}
# 	#normalize each gene
# 	fp_reg = 0.1
# 	#0.1 is defined here because 0.1*mean_num_of_cells_in_cluster
# 	#is epxected to be 3-7, which means that we regulairze
# 	#umicount in the cluster by 3-7.
# 	g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)

# 	return(g_fp_n)
# }