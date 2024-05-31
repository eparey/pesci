import sys
import os
import gc
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
                logger.error(f'Count matrix provided as raw tsv but expected cluster file {clusters} does not exists.')
                sys.exit(1)

        elif ext == '.h5ad':
            fmt = 'h5ad'

        else:
            logger.error(f'{ext} format not supported, please provided .tsv, .h5ad or cellranger directory')
            sys.exit(1)

    elif os.path.isdir(expr_mat):
        fmt = 'cellranger'
        if not os.path.isfile(clusters):
            logger.error(f'Count matrix cellranger directory but expected cluster file {clusters} does not exists.')
            sys.exit(1)

    else:
        logger.error(f'{expr_mat} is not an exisiting file or directory.')
        sys.exit(1)

    return fmt


def load_matrix_tsv(inputfile, cores):
    dt.options.nthreads = cores
    with open(inputfile, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            if len(line) > 2:
                column_ind = line[0].strip('"')
            else:
                logger.error(f'{inputfile} is not tab-delimited.')
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

        logger.info('Count matrix seems to be provided as tab-delimited sparse matrix format, loading with datatable')
        expr = load_matrix_tsv(expr_mat, cores)
        expr = expr[expr.sum(axis = 1) >= 10] #potentially RAN intensive on large datasets, maybe just use scanpy for all

        clust, clust_names = load_cell_clust(clusters,  filter_out_start=filter_out_start)

    elif fmt == 'cellranger':
        logger.info('Count matrix seems to be provided as cellranger directory, loadiing with scanpy')
        sc._settings.ScanpyConfig.n_jobs = cores
        expr = sc.read_10x_mtx(expr_mat)
        sc.pp.filter_genes(expr, min_counts=10, inplace=True) #check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        expr = data.to_df().T

        logger.info('Successfully loaded count matrix with scanpy.')
        clust, clust_names = load_cell_clust(clusters, filter_out_start=filter_out_start)
        logger.info('Successfully loaded clusters')

    elif fmt == 'h5ad':
        logger.info('Count matrix seems to be provided as h5ad, loading with scanpy')
        sc._settings.ScanpyConfig.n_jobs = cores
        data = sc.read(filename = expr_mat)
        if clusters not in data.obs:
            logger.error(f'{clusters} column not found in .obs of h5ad')
        clust = data.obs[clusters].to_dict()
        clust_names = sorted(data.obs[clusters].unique())
        sc.pp.filter_genes(data, min_counts=10, inplace=True) #check it does the same as expr.sum(axis = 1) >= 10, and check we want to do it or not
        expr = data.to_df().T

        logger.info('Successfully loaded count matrix with scanpy.')


    else:
        logger.error(f'{fmt} is not supported.')
        sys.exit(1)

    return expr, clust, clust_names



def load_cell_clust(cells_to_clusters, filter_out_start=None):
    """
    Load cell to clusters
    """
    if not os.path.isfile(cells_to_clusters):
        logger.error(f'{expr_mat} does not exist.')
        sys.exit(1)

    clust = pd.read_csv(cells_to_clusters, sep="\t", header=0, index_col=0)
    clust.dropna(inplace=True)

    if filter_out_start:
        clust = clust[~clust['cluster_name'].astype(str).str.startswith(filter_out_start)]

    val = sorted(clust['cluster_name'].unique())
    clust = clust['cluster_name'].to_dict()
    logger.info(f'Loaded clusters {cells_to_clusters}')
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
        mat = np.array(expr[cells])

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


def normalize_geom_mean_opt(expr, clust, val): #TODO check this gives same as the non-opt

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
        mat = np.array(expr[cells])

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
    data = np.array(data).T + 0.05
    medians = np.median(data, axis=1) + 0.05
    medians = medians.reshape(len(medians), 1)
    data = data / medians
    logger.info(f'{tot} cells, {len(data[:,0])} genes, {len(val)} clusters.')
    return data


def save_outputs(data, genes, clusters, output):
    df = pd.DataFrame(data=data, index=genes, columns=clusters)
    df.to_csv(output, sep='\t')


def normalize_main(expr_mat, cells_to_clusters, output, cores=1, filter_out_start=None):
    fm = validate_input_format(expr_mat, cells_to_clusters)
    matrix, clusters, cluster_names = load_expr_and_clusters(expr_mat, cells_to_clusters, fmt=fm, cores=cores, filter_out_start=filter_out_start)
    norm_matrix = normalize_geom_mean(matrix, clusters, cluster_names)
    save_outputs(norm_matrix, matrix.index, cluster_names, output)
