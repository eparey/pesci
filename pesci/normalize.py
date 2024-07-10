"""
Module with functions to load an RNA-seq single-cell expression count matrix and precomputed
cell clusters and compute normalized fold-change expression per clusters.
"""

import sys
import os
import collections

import logging
import pickle
import gzip
import bz2

import coloredlogs

import datatable as dt
from scipy import sparse
import scanpy as sc
import numpy as np
import pandas as pd

import tqdm

from . import pbar

BAR_FORMAT = pbar.BAR_FORMAT

logger = logging.getLogger(__name__)
coloredlogs.install()


def get_open(ext):

    """
    Infers open function to use for compressed file based on extension.

    Args:
        ext (str): file extention

    Returns:
        function: open function
    """

    if ext in ['.gz', 'gz']:
        return gzip.open
    if ext in ['gz2', '.gz2']:
        return bz2.open

    logger.error('Unsuported format %s', ext)
    sys.exit(1)

def validate_input_format(expr_mat, clusters):

    """
    Infers and validates input file format.
    Infers expression matrix file format (check filename extension, or is directory)
    and check provided cell-to-clusters annotation.

    Three different formats are supported:
        - raw expression matrix as a tab-delimited file
          (.tsv, .txt or .csv, can be compressed in .gz, .gz2)
          + tab-delimited cell-to-clusters file
        - cellranger output directory + tab-delimited cell-to-clusters file
        - scanpy h5ad + name of the column with cluster annotations

    Args:
        expr_mat (str): Expression matrix input file or cellranger directory
        clusters (str): Cell-to-clusters annotation file, or name of the cluster column in the h5ad 

    Returns:
        str: Inferred file format
    """
    open_func = open
    if os.path.isfile(expr_mat):
        name, ext = os.path.splitext(expr_mat)

        if ext in ['.gz', '.gz2']:
            open_func = get_open(ext)
            name, ext = os.path.splitext(name)


        if ext in ['.tsv', '.csv', '.txt']:

            with open_func(expr_mat, 'rt') as infile:
                for line in infile:
                    fmt = 'tsv'
                    if len(line.strip().split(',')) > len(line.strip().split('\t')):
                        fmt = 'csv'
                    break

            if not os.path.isfile(clusters):
                logger.error('Count matrix provided as csv/tsv but expected cluster file %s '
                              'does not exists.', clusters)
                sys.exit(1)

        elif ext == '.h5ad':
            fmt = 'h5ad'

        else:
            logger.error('%s format not supported, please provided a .tsv, .csv, .txt, .h5ad or a '
                         'cellranger directory (.tsv, .csv and .txt can be compressed in .gz '
                         'or .gz2)', ext)
            sys.exit(1)

    elif os.path.isdir(expr_mat):
        fmt = 'cellranger'
        if not os.path.isfile(clusters):
            logger.error('Count matrix cellranger directory but expected cluster file %s '
                         'does not exists.', clusters)
            sys.exit(1)

    else:
        logger.error('%s is not an existing file or directory.', expr_mat)
        sys.exit(1)

    return fmt, open_func


def dt_to_sparse_high_ram(expr, column_ind):
    """
    Converts a datatable Frame count matrix to a scipy sparse matrix. This has high RAM usage since
    the whole matrix is loaded in memory in dense format first.

    Args:
        expr (datatable.Frame): count matrix (with genes still in first column)
        column_inc (str): name of the column with genes

    Returns:
        scipy.sparse.csr_matrix: count matrix in sparse format
    """

    # this loads everything in memory  dense format!
    expr = expr.to_pandas()
    expr.drop(column_ind, axis=1, inplace=True)
    expr = expr.to_numpy().astype(int)
    expr = sparse.csr_matrix(expr)
    return expr

def iterative_dt_to_sparse(dt_expr, cells_per_iter=10000):

    """
    Converts a datatable Frame count matrix to a scipy sparse matrix.
    Conversion is done in chunks of `cells_per_iter` cells (10,000 default) to decrease RAM usage.

    Args:
        expr (datatable.Frame): count matrix (with genes still in first column)
        cells_per_iter (int): size of chunks

    Returns:
        scipy.sparse.csr_matrix: count matrix in sparse format
    """

    logger.info('Converting to a sparse matrix')
    n = len(dt_expr.names)

    #get first 10,000 cells (or all if nb cells<10,000)
    pass_iterative_fill = False
    end = cells_per_iter
    if end > n:
        end = n
        pass_iterative_fill = True

    #dt to dense numpy array
    expr_final = dt_expr[:, 1:end].to_numpy().astype(int)

    #numpy array to sparse
    expr_final = sparse.csr_matrix(expr_final)

    #if matrix has > 10,000 cells we convert to sparse iteratively in chunks of 10,000 cells
    if not pass_iterative_fill:
        expr_subset = [expr_final]
        for i in range(cells_per_iter, n, cells_per_iter):
            end = i + cells_per_iter
            end = min(end, n)
            subset = dt_expr[:, i:end].to_numpy().astype(int)
            subset = sparse.csr_matrix(subset)
            expr_subset.append(subset)
            # expr_final = sparse.hstack((expr_final, subset)) #alternative: hstack in loop
        expr_final = sparse.hstack(expr_subset)
    return expr_final


def load_matrix_tsv(inputfile, cores=1, fmt='tsv', open_func=open):

    """
    Loads expression matrix provided as a text file format, either tab-delimited (fmt='.tsv') or 
    comma delimited (fmt='.csv').

    Args:
        inputfile (str): Name of the input file
        cores (int, optional): Number of cores to use for loading

    Returns:
        tuple: (scipy.sparse array, list, list) expression matrix, genes (rows), cell barcodes (col)
    """
    sep = '\t'
    if fmt == 'csv':
        sep = ','

    dt.options.nthreads = cores
    with open_func(inputfile, 'rt') as infile:
        for line in infile:
            line = line.strip().split(sep)
            if len(line) > 2:
                column_ind = line[0].strip('"')
                if column_ind == '':
                    column_ind = 'C0' #default name given by datatable
            else:
                logger.error('%s does not seem %s delimited.', inputfile, sep)
                sys.exit(1)
            break

    expr = dt.fread(file=inputfile, sep=sep)
    genes = expr[column_ind].to_list()[0]
    cells = list(expr.names[1:])

    #iteratively convert to sparse in chunks of 10,000 cells
    expr = iterative_dt_to_sparse(expr, cells_per_iter=10000)

    return expr, genes, cells


def to_expr_matrix(matrix, genes, cells, clusters):

    """
    Creates a simple named tuple object to represent a cell-by-gene expression matrix.

    Args:
        matrix (numpy.array | scipy.sparse array): expression matrix
        genes (list): list of genes in matrix rows (i.e. rownames)
        cells (lkist): list of cells in matrix columns (i.e. colnames)
        clusters (dict): for each cluster,corresponding cell barcodes (set)

    Returns:
        ExprMatrix namedtuple:
            matrix (scipy.sparse array): Expression matrix
            genes (dict): Names of the genes in rows (for each gene, give its row index)
            cells (dict): Names of the cell in columns (for each cell barcode, give its col index)
            clusters (dict): For each cluster, column indices of corresponding cells (set)
    """

    if len(set(genes)) != len(genes):
        logger.error('Duplicate gene names in count matrix')
        sys.exit(1)

    if len(set(cells)) != len(cells):
        logger.error('Duplicate cell barcodes in count matrix')
        sys.exit(1)

    genes = {gene: idx for idx, gene in enumerate(genes)}
    cells = {cell: idx for idx, cell in enumerate(cells)}
    clusters = {clust: {cells[i] for i in clusters[clust] if i in cells} for clust in clusters}
    clusters = {clust:clusters[clust] for clust in clusters if clusters[clust]} #remove empty clust
    ExprMatrix = collections.namedtuple('ExprMatrix', ['matrix', 'genes', 'cells', 'clusters'])

    return ExprMatrix(matrix, genes, cells, clusters)


def update_dict_of_set(mydict, key, val, filter_out_start=None, keep_only=None):
    """
    Helper function to update a dict of str:set in-place. Used to load cell clusters annotation.
    Skip dict update if cluster name is an empty string or starts with `filter_out_start`.

    Args:
        mydict (dict): dict to update
        key (str): key to update with new value
        val (str): value to add
        filter_out_start (str, optional): ignore clusters whose name start with provided string
    """
    if key == '':
        return
    if filter_out_start and key.startswith(filter_out_start):
        return
    if keep_only and not key.startswith(keep_only):
        return
    mydict[key] = mydict.get(key, set())
    mydict[key].add(val)


def load_expr_and_clusters(expr_mat, clusters,  min_counts=10, fmt='tsv', colclust='cluster_name',
                           cores=1, filter_out_start=None, keep_only=None, broad=None,
                           open_func=open):
    """
    Loads expression matrix and clusters for 3 supported input formats (see validate_input_format).

    Args:
        expr_mat (str): Expression matrix input file or cellranger directory
        clusters (str): Cell-to-clusters annotation file, or name of the cluster column in the h5ad 
        filter_out_start (str, optional): ignore clusters whose name start with provided string
        fmt (str, optional): input format, one of 'tsv' | 'csv' | 'cellranger' | 'h5ad'
        cores (int, optional): number of cores for loading
        colname (str, optional): name of column with clusters in tab-delimited cluster file, default
                                 is 'cluster_name' and will use second column if it does not exist.

    Returns:
        tuple: ExprMat namedtuple stores expression matrix, genes, cells and clusters
               dict with correspondance clusters (key) to broad annoatation (val), empty if
               broad=None
    """

    clust2broad = {}

    #Load tab-delimited count-matrix
    if fmt in ['tsv', 'csv']:

        if fmt == 'tsv':

            logger.info('Count matrix provided as tab-delimited dense matrix format, loading with '
                        'datatable')

        elif fmt == 'csv':
            logger.info('Count matrix provided as comma-delimited dense matrix format, loading with'
                        ' datatable')


        #load matrix
        expr, genes, cells = load_matrix_tsv(expr_mat, cores, fmt, open_func=open_func)

        #filter out genes with no/very low expression
        logger.info('Filtering out lowly-expressed genes (total umi < %s)', min_counts)
        idx_genes_to_keep = np.where(np.sum(expr, axis=1) >= min_counts)[0]
        genes = [genes[i] for i in idx_genes_to_keep]
        expr = expr[idx_genes_to_keep,:]

        #load clusters
        clusters_dict = load_cell_clust(clusters, colname=colclust, keep_only=keep_only,
                                        filter_out_start=filter_out_start)

        if broad:
            clust2broad = load_cell_clust_and_broad(clusters, colname=colclust, colname_broad=broad)

    #Load count matrix in a cellranger output directory
    elif fmt == 'cellranger':
        logger.info('Count matrix provided as cellranger directory, loading with scanpy')
        sc.settings.n_jobs = cores

        #load matrix
        expr = sc.read_10x_mtx(expr_mat)

        #filter out genes with no/very low expression
        logger.info('Filtering out lowly-expressed genes (total umi < %s)', min_counts)
        sc.pp.filter_genes(expr, min_counts=min_counts, inplace=True)
        genes = expr.var.index
        cells = expr.obs.index
        expr = expr.X.T

        logger.info('Successfully loaded count matrix with scanpy.')

        #load clusters
        clusters_dict = load_cell_clust(clusters, colname=colclust,
                                        filter_out_start=filter_out_start, keep_only=keep_only)
        if broad:
            clust2broad = load_cell_clust_and_broad(clusters, colname=colclust, colname_broad=broad)

        logger.info('Successfully loaded clusters')

    #Load count matrix and clusters from .h5ad
    elif fmt == 'h5ad':
        logger.info('Count matrix provided as h5ad, loading with scanpy')
        sc.settings.n_jobs = cores

        #load matrix
        expr = sc.read(filename = expr_mat)

        #TODO check it does the same as expr.sum(axis = 1) >= 10
        sc.pp.filter_genes(expr, min_counts=min_counts, inplace=True)
        logger.info('Filtering out lowly-expressed genes (total umi < %s)', min_counts)

        genes = expr.var.index
        cells = expr.obs.index

        if clusters not in expr.obs:
            logger.error('%s column not found in .obs of h5ad', clusters)
            sys.exit(1)

        #load clusters
        clusters_dict = {}
        for i, row in expr.obs.iterrows():
            update_dict_of_set(clusters_dict, row[clusters], i, filter_out_start, keep_only)

        if broad:
            if broad not in expr.obs:
                logger.error('%s column not found in .obs of h5ad', clusters)
                sys.exit(1)

            clust2broad = {}
            for _, row in expr.obs.iterrows():
                if row[clusters] not in clust2broad:
                    clust2broad[row[clusters]] = row[broad]

        expr = expr.X.T
        logger.info('Successfully loaded count matrix with scanpy.')

    else:
        logger.error('%s is not supported.', fmt)
        sys.exit(1)

    logger.info('%s cells, %s genes, %s clusters.', len(cells), len(genes), len(clusters_dict))

    return to_expr_matrix(expr, genes, cells, clusters_dict), clust2broad


def load_cell_clust(cells_to_clusters, colname='cluster_name', filter_out_start=None,
                    keep_only=None):
    """
    Loads cells cluster annotation from a tab-delimited input file.

    Args:
        cells_to_clusters (str): Input file name
        colname (str, optional): Name of the column with clusters, default is 'cluster_name' and the
                                 second column will be used if 'cluster_name' is not in columns
        filter_out_start (str, optional): ignore clusters whose name start with provided string

    Returns:
        dict: dict of set with key = cluster name, value = set of cell barcodes 
    """
    if not os.path.isfile(cells_to_clusters):
        logger.error('%s does not exist.', cells_to_clusters)
        sys.exit(1)

    clust = {}
    force = False
    sep = '\t'
    with open(cells_to_clusters, 'r', encoding='utf-8') as infile:
        for i, line in enumerate(infile):
            if i == 0:
                line_tmp = [i.strip('"') for i in line.strip().split(sep)]
                if len(line_tmp) == 1:
                    logger.info('Cluster annotation file seems to be comma-delimited %s',
                                cells_to_clusters)
                    sep = ','
                    line = [i.strip('"') for i in line.strip().split(',')]
                else:
                    line = line_tmp

                nb_col_header = len(line)
                if colname not in line:
                    if colname != 'cluster_name':
                        logger.error('Column %s not found in cluster annotation file %s: columms '
                                     'are: %s', colname, cells_to_clusters, ','.join(line))
                        sys.exit(1)
                    else:
                        logger.warning('Using second column of %s as cluster annotation, provide'
                                       ' a column name with --colclust to change this behaviour',
                                        cells_to_clusters)
                        clust_col = 1
                        force = True
                else:
                    clust_col = line.index(colname)
            else:
                line = [i.strip('"') for i in line.strip().split(sep)]
                nb_col = len(line)
                if i == 1:
                    #autodetect if index label is not in header
                    if nb_col == (nb_col_header + 1) and not force:
                        clust_col = clust_col + 1

                if nb_col not in [nb_col_header, nb_col_header + 1]:
                    logger.error('Malformed cluster annotation file %s, '
                                 'different number of columns in header and line %s',
                                 cells_to_clusters, i+1)
                    sys.exit(1)

                cell = line[0]
                clust_curr = line[clust_col]
                update_dict_of_set(clust, clust_curr, cell, filter_out_start, keep_only)

    logger.info('Loaded clusters %s', cells_to_clusters)
    return clust


def load_cell_clust_and_broad(cells_to_clusters, colname='cluster_name',
                              colname_broad='broad_annotation'):
    """
    Loads cluster to broad annotation from a tab-delimited input file.

    Args:
        cells_to_clusters (str): Input file name
        colname_clust (str, optional): Name of the column with clusters, default is 'cluster_name'
                                       and the second column will be used if 'cluster_name' is not
                                       in columns
        colname_broad (str, optional): Name of the column with broad annotation, default is
                                       'broad_annotation' and the third column will be used if
                                       'broad_annotation' is not in columns
    Returns:
        dict: dict of set with key = cluster name, value = broad annotation
    """
    if not os.path.isfile(cells_to_clusters):
        logger.error('%s does not exist.', cells_to_clusters)
        sys.exit(1)

    clust2broad = {}
    force_clust, force_broad = False, False
    sep = '\t'
    with open(cells_to_clusters, 'r', encoding='utf-8') as infile:
        for i, line in enumerate(infile):
            if i == 0:
                line_tmp = [i.strip('"') for i in line.strip().split(sep)]
                if len(line_tmp) == 1:
                    logger.info('Cluster annotation file seems to be comma-delimited %s',
                                cells_to_clusters)
                    sep = ','
                    line = [i.strip('"') for i in line.strip().split(',')]
                else:
                    line = line_tmp

                nb_col_header = len(line)
                if colname not in line:
                    if colname != 'cluster_name':
                        logger.error('Column %s not found in cluster annotation file %s: columms '
                                     'are: %s', colname, cells_to_clusters, ','.join(line))
                        sys.exit(1)
                    else:
                        logger.warning('Using second column of %s as cluster annotation, provide'
                                       ' a column name with --colclust to change this behaviour',
                                        cells_to_clusters)
                        force_clust = True
                else:
                    clust_col = line.index(colname)

                if colname_broad not in line:
                    if colname != 'broad_annotation':
                        logger.error('Column %s not found in cluster annotation file %s: columms '
                                     'are: %s', colname_broad, cells_to_clusters, ','.join(line))
                        sys.exit(1)
                    else:
                        logger.warning('Using third column of %s as broad annotation, provide'
                                       ' a column name with --colbroad to change this behaviour',
                                        cells_to_clusters)
                        broad_col = 2
                        force_broad = True
                else:
                    broad_col = line.index(colname_broad)

            else:
                line = [i.strip('"') for i in line.strip().split(sep)]
                nb_col = len(line)
                if i == 1:
                    #autodetect if index label is not in header
                    if nb_col == (nb_col_header + 1):
                        if not force_clust:
                            clust_col = clust_col + 1
                        if not force_broad:
                            broad_col = broad_col + 1

                if nb_col not in [nb_col_header, nb_col_header + 1]:
                    logger.error('Malformed cluster annotation file %s, '
                                 'different number of columns in header and line %s',
                                 cells_to_clusters, i+1)
                    sys.exit(1)

                clust = line[clust_col]
                if clust not in clust2broad:
                    broad = line[broad_col]
                    clust2broad[clust] = broad
    return clust2broad


def normalize_geom_mean_fc(expr_mat, marker_specificity=0.5, bar_format=None):

    """
    Computes normalized gene expression per cluster

    Args:
        expr_mat (str): Input file name
        marker_specificity (float): marker specificity value, between 0.5 and 0.95
                                    (higher = more specific), used in fold change calculation as
                                    follows: 0.5 computes fold change in cluster a vs median across
                                    clusters, 0.75 computes fold change in cluster a vs 75
                                    expression percentile across clusters (3rd quartile)
        bar_format (str, optional): tqdm bar format, use None for tqdm default

    Returns:
        numpy.array: normalized gene expression per cluster (Fold-change)
    """

    # For each cluster, compute gene expression as geometric mean over cells
    data = []
    all_clusts = sorted(expr_mat.clusters.keys())

    if bar_format:
        bar_format = bar_format.replace('task', 'Computing gene expression per cluster')

    prbar = tqdm.tqdm(all_clusts, colour='#595c79', bar_format=bar_format)
    prbar.unit = ""
    prbar.refresh()
    mean_counts_per_cell = []
    for j, cluster in enumerate(prbar):

        #subset matrix to extract cells of the cluster only
        col_idx_cells = list(expr_mat.clusters[cluster])
        mat = expr_mat.matrix[:,col_idx_cells] #mat is a sparse matrix for RAM efficiency

        #compute mean counts per cell
        mean_counts_per_cell_tmp = np.mean(sparse.csr_matrix.sum(mat, axis=0))
        mean_counts_per_cell.append(mean_counts_per_cell_tmp)

        #compute cluster-level expression per gene using a geometric mean
        mat.data = np.log(mat.data + 1)
        mat = np.array(sparse.csr_matrix.mean(mat, axis=1).T)[0] #todense from here (1 val per gene)
        mat = np.exp(mat) - 1

        #normalize by average counts per cell
        mat = mat / mean_counts_per_cell_tmp

        data.append(mat)

        #update progress bar on success
        if j == len(prbar) - 1:
            prbar.unit = "[done]"
            prbar.refresh()

    # scaling so that 0.05 pseudocount in fold change calc does not completely transform values
    # ideal_count = min(1000, np.median(mean_counts_per_cell))
    data = np.array(data).T * 1000

    # transform expression to fold-change (expression cluster a / median expression across clusters)
    # medians = np.median(data, axis=1)
    medians = np.quantile(data, marker_specificity, axis=1)
    medians = medians.reshape(len(medians), 1)
    data = (data + 0.05) / (medians + 0.05)
    return data


def normalize(expr_mat, cells_to_clusters, output, cores=1, filter_out_start=None,
              keep_only=None, marker_specificity=0.5, colclust='cluster_name', broad=None,
              min_umi=10, bar_format=BAR_FORMAT):
    """
    Funtion to load the data (expression matrix and clusters) and compute per cluster normalized
    gene expression (fold-change).

    Args:
        expr_mat (str): Expression matrix, either a tab-delimited file, cellranger directory or h5ad
        cells_to_clusters (str): input cluster file name or name of cluster column in h5ad
        output (str): output file for nornalized gene-cluster expression matrix
        cores (int, optional): Number of cores to use for loading
        filter_out_start (str, optional): ignore clusters whose name start with provided string
        bar_format (str, optional): tqdm bar format, use None for tqdm default

    """
    fm, open_func = validate_input_format(expr_mat, cells_to_clusters)
    matrix, clust2broad = load_expr_and_clusters(expr_mat, cells_to_clusters, fmt=fm,
                                                 cores=cores, filter_out_start=filter_out_start,
                                                 colclust=colclust, broad=broad,
                                                 open_func=open_func, keep_only=keep_only,
                                                 min_counts=min_umi)
    if broad:
        #this assumes pesci-defined output file name
        outpkl = output.split('_matrix_')[0] + '_clusters_to_broad.pkl'
        with open(outpkl, 'wb') as outfile:
            pickle.dump(clust2broad, outfile)
    norm_matrix = normalize_geom_mean_fc(matrix, marker_specificity=marker_specificity,
                                         bar_format=bar_format)
    df = pd.DataFrame(data=norm_matrix, index=matrix.genes, columns=sorted(matrix.clusters.keys()))
    df.to_csv(output, sep='\t')
