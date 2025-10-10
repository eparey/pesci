"""
Module with functions to load an RNA-seq single-cell expression count matrix and precomputed
cell clusters and compute normalized fold-change expression per cluster.
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

    logger.error('Unsupported format %s', ext)
    sys.exit(1)


def validate_input_format(expr_mat, clusters):

    """
    Infers and validates input file format.
    Infers expression matrix file format (check filename extension, or is_directory)
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
    the whole matrix is loaded in memory in dense format first. Function not used in pesci, retained
    here for testing purposes.

    Args:
        expr (datatable.Frame): count matrix (with genes still in first column)
        column_ind (str): name of the column with genes

    Returns:
        scipy.sparse.csr_matrix: count matrix in sparse format
    """

    # this loads everything in memory  dense format!
    expr = expr.to_pandas()
    expr.drop(column_ind, axis=1, inplace=True)
    expr = expr.to_numpy().astype(int)
    expr = sparse.csr_matrix(expr)
    return expr


def iterative_dt_to_sparse(dt_expr, cells_per_iter=2000):

    """
    Converts a datatable Frame count matrix to a scipy sparse matrix.
    Conversion is done in chunks of `cells_per_iter` cells (2,000 default) to decrease RAM usage.

    Args:
        expr (datatable.Frame): count matrix (with genes still in first column)
        cells_per_iter (int): size of chunks

    Returns:
        scipy.sparse.csr_matrix: count matrix in sparse format
    """

    logger.info('Converting to a sparse matrix')
    n = len(dt_expr.names)

    #get first 10,000 cells (or all if nb cells<2,000)
    pass_iterative_fill = False
    end = cells_per_iter
    if end > n:
        end = n
        pass_iterative_fill = True

    #dt to dense numpy array
    try:
        expr_final = dt_expr[:, 1:end].to_numpy().astype(int)
    except Exception as e:
        if e.__class__.__name__=="TypeError":
            logger.error('The matrix contains non-integer values! Expecting a count matrix, please check!')
            sys.exit(1)

        else:
            raise


    #numpy array to sparse
    expr_final = sparse.csr_matrix(expr_final)

    #if matrix has > 2,000 cells we convert to sparse iteratively in chunks of 2,000 cells
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
        fmt (str, optional): format (tsv or csv)
        open_func (function): open function to use (to handle compressed inputs)

    Returns:
        tuple: (scipy.sparse array, list, list) expression matrix, genes (rows), cell barcodes (col)
    """
    sep = '\t'
    sepname = 'tab'
    if fmt == 'csv':
        sep = ','
        sepname = 'comma'

    dt.options.nthreads = cores
    with open_func(inputfile, 'rt') as infile:
        for line in infile:
            line = line.strip().split(sep)
            lg = len(line)
            if lg > 2:
                column_ind = line[0].strip('"')
                if column_ind == '':
                    column_ind = 'C0' #default name given by datatable
                # colstypes = [dt.str32] + [dt.int32]*(lg-1)
            else:
                logger.error('%s does not seem %s-delimited.', inputfile, sepname)
                sys.exit(1)
            break

    expr = dt.fread(file=inputfile, sep=sep) #columns=colstypes, memory_limit=10000000000
    genes = expr[column_ind].to_list()[0]
    cells = list(expr.names[1:])

    #iteratively convert to sparse in chunks of 2,000 cells
    expr = iterative_dt_to_sparse(expr, cells_per_iter=2000)

    return expr, genes, cells


def to_expr_matrix(matrix, genes, cells, clusters):

    """
    Creates a simple named tuple object to represent a cell-by-gene expression matrix.

    Args:
        matrix (numpy.array | scipy.sparse array): expression matrix
        genes (list): list of genes in matrix rows (i.e. rownames)
        cells (list): list of cells in matrix columns (i.e. colnames)
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


def permute_sparse_matrix(mat, new_row_order):

    """
    Reorders the rows in a scipy sparse matrix using specified index list.

    Args:
        mat (scipy.sparse array): input matrix
        new_row_order (list): new order of row indexes (e.g. [1,0,2,3] swaps first two rows)

    Returns:
        scipy.sparse array: reordered sparse matrix
    """
    new_mat = mat
    i = sparse.eye(mat.shape[0]).tocoo()
    i.row = i.row[new_row_order]
    new_mat = i.dot(new_mat)
    return new_mat


def get_reorder_indexes(prevegenes, newgenes):

    """
    Get correspondence between old gene row indices and desired new ordering.

    Args:
        prevegenes (dict): dict with gene-name:row-index correspondance in old matrix
        newgenes (dict): dict with disered gene-name:row-index correspondance for reoredered matrix

    Returns:
        dict: old_row-indexe:new_row-indexes, sorted by values
    """

    oldgeneidx2new = {prevegenes[gene]:newgenes[gene] for gene in newgenes}
    filtered_out_genesold = {i for i in prevegenes.values() if i not in oldgeneidx2new}
    n = len(oldgeneidx2new)
    j = n
    for i in filtered_out_genesold:
        oldgeneidx2new[i] = j
        j+=1

    oldgeneidx2new = dict(sorted(oldgeneidx2new.items()))
    return oldgeneidx2new, n


def cat_matrices(expr_matrices):

    """
    Concatenate expression matrices (for cases with several sequencing/biological replicates).

    Args:
        expr_matrices (list of ExprMatrix namedtuple): matrices to concatenate

    Returns:
        ExprMatrix namedtuple

    Note:
        see to_expr_matrix() for ExprMatrix object description
    """

    if len(expr_matrices) == 1:
        return expr_matrices[0]

    #update cells indexes
    cellsorder = [cell for i in expr_matrices for cell in i.cells]
    if len(cellsorder) != len(set(cellsorder)):
        logger.error('Duplicate cell barcodes across count matrices')
        sys.exit(1)

    cells = {cell: idx for idx, cell in enumerate(cellsorder)}

    #update clusters to cells dict
    clusters = {}
    for i in expr_matrices:
        #transform previous indexes to new
        oldidx2new = {i.cells[cell]:cells[cell] for cell in i.cells}

        for clust in i.clusters:

            if clust not in clusters:
                clusters[clust] = {oldidx2new[j] for j in i.clusters[clust]}
            else:
                clusters[clust].update({oldidx2new[j] for j in i.clusters[clust]})


    #update genes indexes
    geneorder = sorted(list(set.intersection(*[set(i.genes.keys()) for i in expr_matrices])))
    genes = {gene: idx for idx, gene in enumerate(geneorder)}

    logger.info('Retained %s genes present across all input matrices', len(genes))

    #update matrix by concat all sparse to one (reorder genes first)
    oldgeneidx2new, max_new = get_reorder_indexes(expr_matrices[0].genes, genes)
    matrix = permute_sparse_matrix(expr_matrices[0].matrix, list(oldgeneidx2new.values()))
    mat_to_cat = [matrix[:max_new,:]]
    for i in expr_matrices[1:]:

        oldgeneidx2new, max_new = get_reorder_indexes(i.genes, genes)
        tmp_mat = permute_sparse_matrix(i.matrix, list(oldgeneidx2new.values()))
        mat_to_cat.append(tmp_mat[:max_new,:])

    matrix = sparse.hstack(mat_to_cat)
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
        keep_only (str, optional): keep only cells from clusters starting with provided string
    """

    if key == '':
        return
    if filter_out_start and key.startswith(filter_out_start):
        return
    if keep_only and not key.startswith(keep_only):
        return
    mydict[key] = mydict.get(key, set())
    mydict[key].add(val)


def load_expr_and_clusters(expr_mat, clusters, min_counts=10, fmt='tsv', colclust='cluster_name',
                           cores=1, filter_out_start=None, keep_only=None, broad=None,
                           open_func=open):

    """
    Loads expression matrix and clusters for 3 supported input formats (see validate_input_format).

    Args:
        expr_mat (str): Expression matrix input file or cellranger directory
        clusters (str): Cell-to-clusters annotation file, or name of the cluster column in the h5ad
        min_counts (int, optional): Retain only genes with >= `min_counts` counts
        fmt (str, optional): input format, one of 'tsv' | 'csv' | 'cellranger' | 'h5ad'
        colclust (str, optional): name of column with clusters in tab-delim cluster file, default
                                 is 'cluster_name' and will use second column if it does not exist.
        cores (int, optional): number of cores for loading (default 1)
        filter_out_start (str, optional): ignore clusters whose name start with provided string
        keep_only (str, optional): keep only cells from clusters starting with provided string
        broad (str, optional): name of column with broad clusters in tab-delim cluster file, used
                               to group clusters for plotting, default is to not load broad if not
                               provided.
        open_func (function, optional): auto-detected in pesci, open as file or as compressed file

    Returns:
        tuple:
            ExprMat namedtuple stores expression matrix, genes, cells and clusters
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
        # genes = [lab+'@'+genes[i] for i in idx_genes_to_keep]
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

        # check that values in expr matrix are integers
        subset = expr.X[:10].tocoo().data #look at the 10 first barcodes
        if not all([not (i%1) for i in subset]):
            logger.error('The matrix contains non-integer values! Expecting a count matrix, please check!')
            sys.exit(1)

        # genes = [lab+'@'+g for g in expr.var.index]
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

        #try to find the raw counts in the anndata
        if "counts" in expr.layers: #TODO: could also allow for 'corrected counts' or 'raw'
            expr.X = expr.layers["counts"]
            logger.info('Using count data stored in the "counts" layer.')

        elif expr.raw is not None:
            tmp = expr.raw.to_adata()
            expr.X = tmp[expr.obs_names, expr.var_names].X.copy()
            logger.info('Using count data stored in adata.raw.X.')

        else:
            logger.warning('Pesci needs counts data, but no "counts" layer or adata.raw can be found '
                            'in the h5ad.')
            logger.warning('Using adata.X, but pesci will crash if these are not counts.')


        logger.info('Filtering out lowly-expressed genes (total umi < %s)', min_counts)
        sc.pp.filter_genes(expr, min_counts=min_counts, inplace=True)

        # genes = [lab+'@'+g for g in expr.var.index]
        genes = expr.var.index
        cells = expr.obs.index

        # check that values in expr matrix are integers
        subset = expr.X[:10].tocoo().data #look at the 10 first barcodes
        if not all([not (i%1) for i in subset]):
            logger.error('The matrix contains non-integer values! Expecting a count matrix, please check!')
            sys.exit(1)

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
        keep_only (str, optional): keep only cells from clusters starting with provided string


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


def normalize(expr_mats, cells_to_clusters, output, cores=1, filter_out_start=None,
              keep_only=None, marker_specificity=0.5, colclust='cluster_name', broad=None,
              min_umi=10, bar_format=BAR_FORMAT):

    """
    Funtion to load the data (expression matrix and clusters) and compute per cluster normalized
    gene expression (fold-change).

    Args:
        expr_mat (list or str): Expression matrix (or matrices), either a tab-delimited file,
                                cellranger directory or h5ad
        cells_to_clusters (str): input cluster file name or name of cluster column in h5ad
        output (str): output file for nornalized gene-cluster expression matrix
        cores (int, optional): Number of cores to use for loading
        filter_out_start (str, optional): ignore clusters whose name start with provided string
        keep_only (str, optional): keep only cells from clusters starting with provided string
        marker_specificity (float): marker specificity value, between 0.5 and 0.95
                                    (higher = more specific), used in fold change calculation as
                                    follows: 0.5 computes fold change in cluster a vs median across
                                    clusters, 0.75 computes fold change in cluster a vs 75
                                    expression percentile across clusters (3rd quartile)
        colclust (str, optional): name of column with clusters in tab-delim cluster file, default
                                 is 'cluster_name' and will use second column if it does not exist
        broad (str, optional): name of column with broad clusters in tab-delim cluster file, used
                               to group clusters for plotting, default is to not load broad if not
                               provided.
        min_umi (int, optional): Retain only genes with umi >= `min_umi`
        bar_format (str, optional): tqdm bar format, use None for tqdm default

    """

    all_matrices = []
    if type(expr_mats) == str:
        expr_mats = [expr_mats]
    for expr_mat in expr_mats:
        fm, open_func = validate_input_format(expr_mat, cells_to_clusters)
        matrix, clust2broad = load_expr_and_clusters(expr_mat, cells_to_clusters, fmt=fm,
                                                     cores=cores, filter_out_start=filter_out_start,
                                                     colclust=colclust, broad=broad,
                                                     open_func=open_func, keep_only=keep_only,
                                                     min_counts=min_umi)
        all_matrices.append(matrix)

    matrix = cat_matrices(all_matrices)

    if broad:
        #FIXME this assumes pesci-defined output file name
        outpkl = output.split('_matrix_')[0] + '_clusters_to_broad.pkl'
        with open(outpkl, 'wb') as outfile:
            pickle.dump(clust2broad, outfile)
    norm_matrix = normalize_geom_mean_fc(matrix, marker_specificity=marker_specificity,
                                         bar_format=bar_format)
    df = pd.DataFrame(data=norm_matrix, index=matrix.genes, columns=sorted(matrix.clusters.keys()))
    df.to_csv(output, sep='\t')

