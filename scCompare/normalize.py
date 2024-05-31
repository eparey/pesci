import logging
import coloredlogs

import datatable as dt
import numpy as np
import pandas as pd

import tqdm


logger = logging.getLogger(__name__)
coloredlogs.install()


def load_gene_cell_expr(expr_mat, cores=1):
    """
    Load gene-cell expression matrix (cells are columns, rows genes).
    """

    dt.options.nthreads = cores

    with open(expr_mat, 'r') as infile:
        for line in infile:
            column_ind = line.strip().split('\t')[0].strip('"')
            break

    expr = dt.fread(file=expr_mat, sep='\t')
    # expr = pd.read_csv(expr_mat, sep="\t", engine='pyarrow', header=0, index_col=0)
    expr = expr.to_pandas()
    expr.set_index(column_ind, inplace=True)

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
    clust.dropna(inplace=True)
    if filter_out_start:
        clust = clust[~clust['cluster_name'].astype(str).str.startswith(filter_out_start)]
    clust = clust['cluster_name'].to_dict()
    val = sorted({str(i) for i in clust.values()})
    return clust, val



def geometric_mean(vector):
    return (np.exp(np.mean(np.log(1 + vector)))) - 1


def normalize_geom_mean(expr, clust, val):

    # For each cluster, compute gene expression as geometric mean over cells (do parallel?)
    data = []
    tot = 0
    pbar = tqdm.tqdm(val, colour='#595c79', bar_format='{percentage:3.0f}% |{bar:50}| Computing expression per cluster \x1B[1;32m{unit}')
    pbar.unit = ""
    pbar.refresh()
    for j, v in enumerate(pbar):
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

        if j == len(pbar) - 1:
            pbar.unit = "[done]"
            pbar.refresh()

    # Normalize gene expression by dividing by its median of expression across clusters
    data = np.array(data).T
    medians = np.median(data, axis=1)
    for i in range(len(medians)):
        data[i,:] = (data[i,:] + 0.05) / (medians[i] + 0.05) #if we keep the code, remove the loop and do matrix operation

    logger.info(f'{tot} cells, {len(data[:,0])} genes, {len(val)} clusters.')
    return data


def save_outputs(data, genes, clusters, output):
    df = pd.DataFrame(data=data, index=genes, columns=clusters)
    df.to_csv(output, sep='\t')


def normalize_main(expr_mat, cells_to_clusters, output, cores=1, filter_out_start=None):
    matrix = load_gene_cell_expr(expr_mat, cores=cores)
    clusters, cluster_names = load_cell_clust(cells_to_clusters, filter_out_start=filter_out_start)
    norm_matrix = normalize_geom_mean(matrix, clusters, cluster_names)
    save_outputs(norm_matrix, matrix.index, cluster_names, output)
