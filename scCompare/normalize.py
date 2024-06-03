import sys
import os
import logging
import coloredlogs

import datatable as dt
import scanpy as sc
import numpy as np
import pandas as pd

import tqdm

logger = logging.getLogger(__name__)
coloredlogs.install()

def validate_input_format(expr_mat, clusters):

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


def load_matrix_tsv(inputfile, cores):
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
    # expr = pd.read_csv(expr_mat, sep="\t", engine='pyarrow', header=0, index_col=0)
    expr = expr.to_pandas()
    expr.set_index(column_ind, inplace=True)
    return expr


def load_expr_and_clusters(expr_mat, clusters, filter_out_start=None, fmt='tsv', cores=1):
    """
    Load gene-cell expression matrix (uses R convention where cells are columns, rows are genes).


    #Filter out non-expressed genes
    # expr = expr.loc[~(expr==0).all(axis=1)]
    """

    if fmt == 'tsv':

        logger.info('Count matrix provided as tab-delimited dense matrix format, loading with '
                    'datatable')
        expr = load_matrix_tsv(expr_mat, cores)
        #TODO potentially RAN intensive on large datasets, maybe just use scanpy for all
        expr = expr[expr.sum(axis = 1) >= 10]

        clust, clust_names = load_cell_clust(clusters,  filter_out_start=filter_out_start)

    elif fmt == 'cellranger':
        logger.info('Count matrix provided as cellranger directory, loadiing with scanpy')
        sc.settings.n_jobs = cores
        expr = sc.read_10x_mtx(expr_mat)

        #TODO check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        sc.pp.filter_genes(expr, min_counts=10, inplace=True)
        expr = expr.to_df().T

        logger.info('Successfully loaded count matrix with scanpy.')
        clust, clust_names = load_cell_clust(clusters, filter_out_start=filter_out_start)
        logger.info('Successfully loaded clusters')

    elif fmt == 'h5ad':
        logger.info('Count matrix provided as h5ad, loading with scanpy')
        sc.settings.n_jobs = cores
        expr = sc.read(filename = expr_mat)
        if clusters not in expr.obs:
            logger.error('%s column not found in .obs of h5ad', clusters)
        clust = expr.obs[clusters].to_dict()
        clust_names = sorted(expr.obs[clusters].unique())

        #TODO check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        sc.pp.filter_genes(expr, min_counts=10, inplace=True)
        expr = expr.to_df().T

        logger.info('Successfully loaded count matrix with scanpy.')


    else:
        logger.error('%s is not supported.', fmt)
        sys.exit(1)

    return expr, clust, clust_names



def load_cell_clust(cells_to_clusters, filter_out_start=None):
    """
    Load cell to clusters
    """
    if not os.path.isfile(cells_to_clusters):
        logger.error('%s does not exist.', cells_to_clusters)
        sys.exit(1)

    clust = pd.read_csv(cells_to_clusters, sep="\t", header=0, index_col=0)
    clust.dropna(inplace=True)

    if filter_out_start:
        clust = clust[~clust['cluster_name'].astype(str).str.startswith(filter_out_start)]

    val = sorted(clust['cluster_name'].unique())
    clust = clust['cluster_name'].to_dict()
    logger.info('Loaded clusters %s', cells_to_clusters)
    return clust, val


def normalize_geom_mean_fc(expr, clust, val): #TODO check this gives same as the non-opt

    # For each cluster, compute gene expression as geometric mean over cells (do parallel?)
    data = []
    tot = 0
    #TODO global def of pbar format + check color in  non-colored term
    pbar = tqdm.tqdm(val, colour='#595c79', bar_format='{percentage:3.0f}% |{bar:50}| Computing expression per cluster \x1B[1;32m{unit}')
    pbar.unit = ""
    pbar.refresh()
    mean_sizes = []
    for j, v in enumerate(pbar):
        cells = [i for i in clust if clust[i]==v]
        lg = len(cells)
        tot += lg
        col_idx = expr.columns.get_indexer(cells)
        mat = np.array(expr.values[:,col_idx])

        mean_clustsize = np.mean(np.sum(mat, axis=0))
        mean_sizes.append(mean_clustsize)

        mat = np.log(mat + 1)
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
    logger.info('%s cells, %s genes, %s clusters.', tot, len(data[:,0]), len(val))
    return data


def save_outputs(data, genes, clusters, output):
    df = pd.DataFrame(data=data, index=genes, columns=clusters)
    df.to_csv(output, sep='\t')


def normalize_main(expr_mat, cells_to_clusters, output, cores=1, filter_out_start=None):
    fm = validate_input_format(expr_mat, cells_to_clusters)
    matrix, clusters, cluster_names = load_expr_and_clusters(expr_mat, cells_to_clusters, fmt=fm,
                                                             cores=cores,
                                                             filter_out_start=filter_out_start)
    norm_matrix = normalize_geom_mean_fc(matrix, clusters, cluster_names)
    save_outputs(norm_matrix, matrix.index, cluster_names, output)
