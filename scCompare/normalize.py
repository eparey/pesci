import sys
import os
import collections

import logging
import coloredlogs

import datatable as dt
from scipy import sparse
import scanpy as sc
import numpy as np
import pandas as pd

import tqdm

logger = logging.getLogger(__name__)
coloredlogs.install()


def to_expr_matrix(matrix, genes, cells, clusters):

    """
    Simple named tuple object to represent a cell-by-gene expression matrix.

    Args:


    Returns:
        ExprMatrix namedtuple:
            matrix (scipy.sparse array): Expression matrix
            genes (dict): Names of the genes in rows (for each gene, give its row index)
            cells (dict): Names of the cell in columns (for each cell barcode, give its col index)
            clusters (dict): For each cluster, column indices of corresponding cells (set)
    """

    genes = {gene: idx for idx, gene in enumerate(genes)}
    cells = {cell: idx for idx, cell in enumerate(cells)}
    clusters = {clust: {cells[i] for i in clusters[clust]} for clust in clusters}

    ExprMatrix = collections.namedtuple('ExprMatrix', ['matrix', 'genes', 'cells', 'clusters'])

    return ExprMatrix(matrix, genes, cells, clusters)


def validate_input_format(expr_mat, clusters):

    """
    Infers and validates input file format.
    Infers expression matrix file format (check filename extension, or is directory)
    and check provided cell-to-clusters annotation.

    Three different formats are supported:
        - raw expression matrix as a tab-delimited file
          (.tsv, .txt or .csv, can be compressed in .gz, .zip, .gz2, .tar, .tgz)
          + tab-delimited cell-to-clusters file
        - cellranger output directory + tab-delimited cell-to-clusters file
        - scanpy h5ad + name of the column with cluster annotations

    Args:
        expr_mat (str): Expression matrix input file or cellranger directory
        clusters (str): Cell-to-clusters annotation file, or name of the cluster column in the h5ad 

    Returns:
        str: Inferred file format
    """

    if os.path.isfile(expr_mat):
        name, ext = os.path.splitext(expr_mat)

        if ext in ['.gz', '.zip', '.gz2', '.tar', '.tgz']:
            name, ext = os.path.splitext(name)

        if ext in ['.tsv', '.csv', '.txt']:
            fmt = 'tsv'

            if not os.path.isfile(clusters):
                logger.error('Count matrix provided as tsv but expected cluster file %s '
                              'does not exists.', clusters)
                sys.exit(1)

        elif ext == '.h5ad':
            fmt = 'h5ad'

        else:
            logger.error('%s format not supported, please provided .tsv, .h5ad or cellranger '
                         'directory', ext)
            sys.exit(1)

    elif os.path.isdir(expr_mat):
        fmt = 'cellranger'
        if not os.path.isfile(clusters):
            logger.error('Count matrix cellranger directory but expected cluster file %s '
                         'does not exists.', clusters)
            sys.exit(1)

    else:
        logger.error('%s is not an exisiting file or directory.', expr_mat)
        sys.exit(1)

    return fmt


def load_matrix_tsv(inputfile, cores=1):

    """
    Loads expression matrix provided as a tab-delimited file using datatable.

    Args:
        inputfile (str): Name of the input file
        cores (int, optional): Number of cores to use for loading

    Returns:
        (scipy.sparse array, list, list): expression matrix, gene names (rows), cell barcodes (col)
    """

    dt.options.nthreads = cores
    with open(inputfile, 'r', encoding = "utf-8") as infile:
        for line in infile:
            line = line.strip().split('\t')
            if len(line) > 2:
                column_ind = line[0].strip('"')
            else:
                logger.error('%s is not tab-delimited.', inputfile)
                sys.exit(1)
            break

    expr = dt.fread(file=inputfile, sep='\t')
    expr = expr.to_pandas()
    genes = expr[column_ind].tolist()
    expr.drop(column_ind, axis=1, inplace=True)
    cells = list(expr)
    expr = expr.to_numpy()
    expr = sparse.csr_matrix(expr)
    return expr, genes, cells


def load_expr_and_clusters(expr_mat, clusters, filter_out_start=None, fmt='tsv', cores=1):
    """
    Loads expression matrix and clusters for 3 supported input formats (see validate_input_format).

    Args:
        expr_mat (str): Expression matrix input file or cellranger directory
        clusters (str): Cell-to-clusters annotation file, or name of the cluster column in the h5ad 
        filter_out_start (str, optional): ignore clusters whose name start with provided string
        fmt (str, optional): input format, one of 'tsv' | 'cellranger' | 'h5ad'
        cores (int, optional): number of cores for loading

    Returns:
        (ExprMat): object storing expression matrix, genes, cells and clusters
    """

    if fmt == 'tsv':

        logger.info('Count matrix provided as tab-delimited dense matrix format, loading with '
                    'datatable')

        expr, genes, cells = load_matrix_tsv(expr_mat, cores)
        idx_genes_to_keep = np.where(np.sum(expr, axis=1) >= 10)[0]

        genes = [genes[i] for i in idx_genes_to_keep]
        expr = expr[idx_genes_to_keep,:]
        clusters_dict = load_cell_clust(clusters,  filter_out_start=filter_out_start)

    elif fmt == 'cellranger':
        logger.info('Count matrix provided as cellranger directory, loading with scanpy')
        sc.settings.n_jobs = cores
        expr = sc.read_10x_mtx(expr_mat)

        #TODO check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        sc.pp.filter_genes(expr, min_counts=10, inplace=True)
        genes = expr.var
        cells = expr.obs
        expr = expr.X.T

        logger.info('Successfully loaded count matrix with scanpy.')
        clusters_dict = load_cell_clust(clusters, filter_out_start=filter_out_start)
        logger.info('Successfully loaded clusters')

    elif fmt == 'h5ad':
        logger.info('Count matrix provided as h5ad, loading with scanpy')
        sc.settings.n_jobs = cores
        expr = sc.read(filename = expr_mat)

        #TODO check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        sc.pp.filter_genes(expr, min_counts=10, inplace=True)
        genes = expr.var.index
        cells = expr.obs.index

        if clusters not in expr.obs:
            logger.error('%s column not found in .obs of h5ad', clusters)

        clusters_dict = {}
        for i, row in expr.obs.iterrows():
            update_dict_of_set(clusters_dict, row[clusters], i, filter_out_start)

        expr = expr.X.T
        logger.info('Successfully loaded count matrix with scanpy.')

    else:
        logger.error('%s is not supported.', fmt)
        sys.exit(1)
    return to_expr_matrix(expr, genes, cells, clusters_dict)


def update_dict_of_set(mydict, key, val, filter_out_start=None):
    if key == '':
        return
    if filter_out_start and key.startswith(filter_out_start):
        return
    mydict[key] = mydict.get(key, set())
    mydict[key].add(val)


def load_cell_clust(cells_to_clusters, colname='cluster_name', filter_out_start=None):
    """
    Load cell to clusters
    """
    if not os.path.isfile(cells_to_clusters):
        logger.error('%s does not exist.', cells_to_clusters)
        sys.exit(1)

    clust = {}
    with open(cells_to_clusters, 'r', encoding='utf-8') as infile:
        for i, line in enumerate(infile):
            line = [i.strip('"') for i in line.strip().split('\t')]
            if i == 0:
                if colname not in line:
                    logger.error('Column %s not found in cluster annotation file %s', colname,
                                 cells_to_clusters)
                    sys.exit(1)
                clust_col = line.index(colname)
            else:
                cell = line[0]
                clust_curr = line[clust_col]
                update_dict_of_set(clust, clust_curr, cell, filter_out_start)

    logger.info('Loaded clusters %s', cells_to_clusters)
    return clust


def normalize_geom_mean_fc(expr_mat):

    # For each cluster, compute gene expression as geometric mean over cells (do parallel?)
    data = []
    tot = 0

    all_clusts = sorted(expr_mat.clusters.keys())

    #TODO global def of pbar format + check color in  non-colored term
    pbar = tqdm.tqdm(all_clusts, colour='#595c79', bar_format='{percentage:3.0f}% |{bar:50}| Computing expression per cluster \x1B[1;32m{unit}')
    pbar.unit = ""
    pbar.refresh()
    mean_sizes = []
    for j, cluster in enumerate(pbar):
        col_idx_cells = list(expr_mat.clusters[cluster])
        tot += len(col_idx_cells)
        mat = expr_mat.matrix[:,col_idx_cells]

        mean_clustsize = np.mean(sparse.csr_matrix.sum(mat, axis=0))
        mean_sizes.append(mean_clustsize)

        mat = np.log(mat.toarray() + 1)
        mat = np.exp(np.mean(mat, axis=1)) - 1
        mat = mat / mean_clustsize

        data.append(mat)

        if j == len(pbar) - 1:
            pbar.unit = "[done]"
            pbar.refresh()

    # Normalize gene expression by dividing by its median of expression across clusters
    ideal_cell_size = min(1000, np.median(mean_sizes))
    data = np.array(data).T * ideal_cell_size
    medians = np.median(data, axis=1)
    medians = medians.reshape(len(medians), 1)
    data = (data + 0.05) / (medians + 0.05)
    logger.info('%s cells, %s genes, %s clusters.', tot, len(data[:,0]), len(expr_mat.clusters))
    return data


def normalize_main(expr_mat, cells_to_clusters, output, cores=1, filter_out_start=None):
    fm = validate_input_format(expr_mat, cells_to_clusters)
    matrix = load_expr_and_clusters(expr_mat, cells_to_clusters, fmt=fm,
                                                             cores=cores,
                                                             filter_out_start=filter_out_start)
    norm_matrix = normalize_geom_mean_fc(matrix)
    df = pd.DataFrame(data=norm_matrix, index=matrix.genes, columns=sorted(matrix.clusters.keys()))
    df.to_csv(output, sep='\t')
