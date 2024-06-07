"""
Module with functions to compute expression conservation scores for 1-to-1 orthologs across
scRNA-seq datasets in two species, using the ICC procedure. Also selects best gene pairs
for many-to-many or 1-to-many orthologs (pair with the highest co-expression conservation).
"""

import sys
import logging

import itertools

import multiprocessing
import traceback
import signal

import coloredlogs

import tqdm

import numpy as np
import pandas as pd

import networkx as nx

from . import pbar
from . import normalize as nm

BAR_FORMAT = pbar.BAR_FORMAT

logger = logging.getLogger(__name__)
coloredlogs.install()


def init_ec(mat1, mat2):

    """
    Initializes the expression conservation score (EC score), using unweighted orthologs
    co-correlation.

    mat1 and mat2 should be expression correlation matrices of the same size (ngenes * ngenes) with
    genes in the same order in both matrices.

    Args:
        mat1, mat2 (numpy.array): orthologs correlation matrix from species 1 and 2 respectively

    Returns:
        numpy.array: vector with EC scores
    """

    ec = []
    n = len(mat1[:,0])
    for i in range(n):
        ec.append(max(np.corrcoef(mat1[i,:], mat2[i,:])[0,1], 0)) #co-correlation rows i or 0 if neg
    return np.array(ec)


def init_ec_optimized_einsum(mat1, mat2):

    """
    Initializes the expression conservation score (EC score), using unweighted orthologs
    co-correlation. This is the same as init_ec() but using numpy einsum for better performance
    (to the cost of a decreased readability).

    mat1 and mat2 should be expression correlation matrices of the same size (ngenes * ngenes) with
    genes in the same order in both matrices.

    Args:
        mat1, mat2 (numpy.array): orthologs correlation matrix from species 1 and 2 respectively

    Returns:
        numpy.array: vector with EC scores
    """

    n = len(mat1[:,0])

    #center matrix (mat - mat average along rows)
    mat1_m = mat1 - np.einsum("ij->i", mat1, optimize='optimal') / np.double(n)
    mat2_m = mat2 - np.einsum("ij->i", mat2, optimize='optimal') / np.double(n)

    #get the diagonal of the covariance matrix (as we only need the co-correlation of orthologs)
    cov = np.einsum("ij,ij->j", mat1_m, mat2_m, optimize='optimal')

    #get the variance
    var1 = np.einsum("ij,ij->j", mat1_m, mat1_m, optimize='optimal')
    var2 = np.einsum("ij,ij->j", mat2_m, mat2_m, optimize='optimal')

    #get variance product
    varprod = np.einsum("i,i->i", var1, var2, optimize='optimal')

    #correlation coefficient
    ec = cov / np.sqrt(varprod)

    return ec.clip(min=0) #transform negative corr to 0


def update_ec(mat1, mat2, prev_ec):

    """
    Updates the expression conservation score (EC score), using the ICC procedure
    (Iterative Correlation of coexpression).

    mat1 and mat2 should be expression correlation matrices of the same size (ngenes * ngenes) with
    genes in the same order in both matrices (and same order as the ec vector).

    Args:
        mat1, mat2 (numpy.array): orthologs correlation matrix from species 1 and 2 respectively
        prev_ec (numpy.array): ec score vector from previous iteration

    Returns:
        numpy.array: vector with EC scores
    """

    def __wcov(v1, v2, w):
        r = np.sum(w * (v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w)))
        return r / np.sum(w)
        #np.average((v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w)), weights=w)

    def __wcorr(v1, v2, w):
        res = __wcov(v1, v2, w) / np.sqrt(__wcov(v1, v1, w) * __wcov(v2, v2, w))
        return res

    n = len(mat1[:,0])
    ec = []
    for i in range(n):
        ec.append(max(__wcorr(mat1[i,:], mat2[i,:], prev_ec), 0))
    return np.array(ec)


def update_ec_optimized_einsum(mat1, mat2, prev_ec):

    """
    Update the expression conservation score (EC score) from a previous iteration.
    Same as update_ec() but using numpy einsum for better performance (to the cost of a decreased
    readability).

    mat1 and mat2 should be expression correlation matrices of the same size (ngenes * ngenes) with
    genes in the same order in both matrices (and same order as the ec vector).

    Args:
        mat1, mat2 (numpy.array): orthologs correlation matrix from species 1 and 2 respectively
        prev_ec (numpy.array): ec score vector from previous iteration

    Returns:
        numpy.array: vector with EC scores
    """

    n = len(mat1[:,0])

    #repeat weight vector n times for matrix operations
    wmat = np.tile(prev_ec, (n, 1))

    #weighted centering of matrices
    mat1_m = mat1 - np.einsum("ji,ij->j", wmat, mat1, optimize='optimal') / sum(prev_ec)
    mat2_m = mat2 - np.einsum("ji,ij->j", wmat, mat2, optimize='optimal') / sum(prev_ec)

    #weighted covariance matrix
    #get only diagonal of the covariance matrix (as we only need the co-correlation of orthologs)
    cov12 = np.einsum("ij,ij->ij", mat1_m, mat2_m, optimize='optimal')
    wcov12 = np.einsum("ji,ij->j", wmat, cov12, optimize='optimal') / sum(prev_ec)

    cov1 = np.einsum("ij,ij->ij", mat1_m, mat1_m, optimize='optimal')
    wcov1 = np.einsum("ji,ij->j", wmat, cov1, optimize='optimal') / sum(prev_ec)

    cov2 = np.einsum("ij,ij->ij", mat2_m, mat2_m, optimize='optimal')
    wcov2 = np.einsum("ji,ij->j", wmat, cov2, optimize='optimal') / sum(prev_ec)

    prod = np.einsum("i,i->i", wcov1, wcov2, optimize='optimal')

    #weighted Pearson correlation coefficent
    res = wcov12 / np.sqrt(prod)
    res[np.isnan(res)] = 0
    return res.clip(min=0) #transform negative corr to 0


def update_ec_until_convergence(mat1, mat2, ec, max_iter=100):

    """
    Iteratively updates the expression conservation score (EC score) using the ICC procedure
    (Iterative Correlation of coexpression) until convergence or max_iter is reach.

    mat1 and mat2 should be expression correlation matrices of the same size (ngenes * ngenes) with
    genes in the same order in both matrices (and same order as the ec vector).

    Args:
        mat1, mat2 (numpy.array): orthologs correlation matrix from species 1 and 2 respectively
        ec (numpy.array): initialized ec score vector
        max_iter (int, optional): maximum iteration

    Returns:
        numpy.array: vector with EC scores
    """

    i = 0
    diff = np.Inf
    while i < max_iter and diff > 0.05:
        ec_new = update_ec_optimized_einsum(mat1, mat2, ec)
        diff = np.sum((ec - ec_new)**2)
        ec = ec_new
        i += 1

    if i == max_iter:
        logger.warning('ICC convergence not reached, consider increasing max_iter.')

    return ec


def compute_expression_conservation(matrix_sp_a, matrix_sp_b):

    """
    Compute expression conservation score (EC score) using the ICC procedure
    (Iterative Correlation of coexpression) until convergence or max_iter is reach.

    mat1 and mat2 should be expression matrices with the same number of rows (n genes) with
    genes in the same order in both matrices.

    Args:
        matrix_sp_a, matrix_sp_b (numpy.array): expression matrix from species 1 and 2 respectively

    Returns:
        numpy.array: vector with EC scores
    """

    #transform each expression matrix into a correlation co-expression matrix
    co_a = np.corrcoef(matrix_sp_a)
    co_b = np.corrcoef(matrix_sp_b)

    #for each row in a, compute correlation with the corresponding row in b to init EC
    expr_conservation = init_ec_optimized_einsum(co_a, co_b)

    #Iteratively update EC using weighted correlation until convergence or max iter is reached
    expr_conservation_final = update_ec_until_convergence(co_a, co_b, expr_conservation)

    return expr_conservation_final


def column_wise_wcorr_einsum(mat1, mat2, w):

    """
    Computes weighted correlation for all pairs of columns mat1 - mat2, using numpy einsum for
    performance.

    mat1 and mat2 are expected to be expression matrix of the form (n genes x k clusters) and
    (n genes x l clusters), respectively, with genes in the same order (= same order in weight
    vector w).

    Args:
        mat1, mat2 (numpy.array): orthologs expression matrix from species 1 and 2 respectively
        w (numpy.array): weight vector

    Returns:
        numpy.array: matrix (k, l) weighted pearson correlation coefficient for all column pairs
    """

    #n variable is unused but it helps readability
    (n, k) = mat1.shape  # n genes k clust
    (n, l) = mat2.shape  # n genes l clust

    wmat1 = np.tile(w, (k, 1)) #repeat weight vector k times for matrix operations
    wmat2 = np.tile(w, (l, 1)) #repeat weight vector l times for matrix operations

    #weighted centering of matrix
    mat1_m = mat1 - np.einsum("kn,nk->k", wmat1, mat1, optimize='optimal') / sum(w)
    mat2_m = mat2 - np.einsum("ln,nl->l", wmat2, mat2, optimize='optimal') / sum(w)

    #weighted covariance matrix
    wmat1_m = np.einsum("nk,kn->nk", mat1_m, wmat1, optimize='optimal')
    tmp12 = np.einsum("nk,nl->kl", wmat1_m, mat2_m, optimize='optimal')
    wcov12 = tmp12 / sum(w)

    #weighted variance
    cov1 = np.einsum("nk,nk->nk", mat1_m, mat1_m, optimize='optimal')
    wcov1 = np.einsum("kn,nk->k", wmat1, cov1, optimize='optimal') / sum(w)

    #weighted variance
    cov2 = np.einsum("nl,nl->nl", mat2_m, mat2_m, optimize='optimal')
    wcov2 = np.einsum("ln,nl->l", wmat2, cov2, optimize='optimal') / sum(w)

    #product of variances
    prod = np.einsum("k,l->kl", wcov1, wcov2, optimize='optimal')

    #Pearson weighted correlation coefficent
    res = wcov12 / np.sqrt(prod)

    return res.clip(min=0)


def load_orthologs(input_file, genes_sp_a, genes_sp_b):

    """
    Load orthologs gene pairs from tab-delimited input file (one pair per line). 

    Args:
        mat1, mat2 (numpy.array): orthologs expression matrix from species 1 and 2 respectively
        w (numpy.array): weight vector

    Returns:
        (list, list): list of one2one orthologs and list many-to-many / many-to-one orthologs, 
                      each list is a list of tuple (genes in sp1, genes in species 2).
    """

    #load all orthologous pairs
    edges = []
    with open(input_file, 'r', encoding = "utf-8") as infile:
        for line in infile:
            g1, g2 = line.strip().split('\t')
            if g1 in genes_sp_a and g2 in genes_sp_b:
                edges.append((g1, g2))
            elif g2 in genes_sp_a and g1 in genes_sp_b:
                edges.append((g2, g1))

    #transform list of orthologous pairs to graph
    homologs_graph = nx.Graph()
    homologs_graph.add_edges_from(edges)

    #one-to-one orthologs will be connected to each other only --> connected component of the graph
    components = nx.connected_components(homologs_graph)
    one2one, many_ortho = [], []

    #go through connected component to extract one-to-on and many-to-many / many-to-one
    for component in components:
        genes = set(component)
        if len(component) == 2:
            g1, g2 = genes
            if g1 in genes_sp_b:
                g1, g2 = g2, g1
            one2one.append((g1, g2))

        else:
            tmp_genes_a = genes.intersection(genes_sp_a)
            tmp_genes_b = genes.intersection(genes_sp_b)
            many_ortho.append((tuple(tmp_genes_a), tuple(tmp_genes_b)))

    return one2one, many_ortho


def load_cluster_expression_matrix(input_file):

    """
    Load tab-delimited cluster gene expression matrix.

    Args:
        inputfile (str): name of the expression matrix file

    Returns:
        ExprMat namedtuple: stores expression matrix, genes and cell_types (clusters)
    """

    #load data
    df = pd.read_csv(input_file, sep="\t", header=0, index_col=0)

    #remove constant rows
    unique_values = df.nunique(axis=1)
    constant_rows = unique_values[unique_values == 1].index.tolist()
    df = df.drop(constant_rows)

    #get gene names
    genes = df.index.values.tolist()

    #get cell types
    cell_types = df.columns.tolist()


    expr = nm.to_expr_matrix(np.array(df), genes, cell_types, {})

    return expr


def init_worker():
    """
    Initialisation of worker for parallel operations
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def worker_add_gene_pairs(manyortho, a, b, mat1, mat2, max_combin):

    """
    Selects the best gene pairs (most conserved expression) for many-to-many or 1-to-many orthologs.

    Args:
        manyortho (list of tuples): lits of sp1 - sp2 gene pairs that are many to many orthologs.
        a, b (numpy.array): gene - cluster expression matrix for the randomly selected subset
                            of 1-to-1 orthologs, for sp1 and sp2 respectively.
        mat1, mat2 (tuple): tuple with gene - cluster expression matrix + list of genes (rows)
                            for sp1 and sp2, respectively
        max_combin (int): maximum accepted number or pairwise combinations, massively multigeneic
                          families will be skipped.

    Returns:
        list of str: list of size two, (1) tab-separated tested gene pairs (2) tab-seprated 
                     corresponding scores 
    """

    try:

        res = []
        mat1, genes1 = mat1
        mat2, genes2 = mat2
        for og in manyortho:

            #list all orthologs pairings (create list of tuples)
            combin = [(i, j) for (i, j) in list(itertools.product(og[0], og[1])) if i in genes1
                      and j in genes2]

            if len(combin) > max_combin:
                res.append(f'skipped {og} - {len(og[0])} sp1 genes - {len(og[1])} sp2 genes - '
                           f'{len(combin)} combinations\n')

            else:
                #add all sp1 genes to matrix1 (with duplicates to accomodate all possible pairings)
                manyortho_a = [genes1[i[0]] for i in combin]

                a_w_manyortho_tmp = np.concatenate((a, mat1[manyortho_a,:]), axis=0)

                #add all sp2 genes to matrix2 (with duplicates to accomodate all possible pairings)
                manyortho_b = [genes2[i[1]] for i in combin]

                b_w_manyortho_tmp = np.concatenate((b, mat2[manyortho_b,:]), axis=0)

                #compute expression conservation for all possible pairings
                expr_cons_w_manyortho_current = compute_expression_conservation(a_w_manyortho_tmp,
                                                                                b_w_manyortho_tmp)

                scores = '\t'.join([str(i) for i in expr_cons_w_manyortho_current[1000:]])
                tmp = '\t'.join(['+'.join(comb) for comb in combin])+ '\n' + scores +'\n'
                res.append(tmp)

        return res

    except Exception:
        traceback.print_exc()
        raise


def write_ec_one2one(one2one, expr_conservation_orthologs, out_ortho):

    """
    Writes scores to a tab-seprated file.

    Args:
        ono2one (list): gene names (column names of expr_conservation_orthologs)
        expr_conservation_orthologs (numpy.array): vector with EC scores for each gene, genes in
                                                   the same order as in one2one
        out_ortho (int): name of output file
    """

    with open(out_ortho, 'w', encoding = "utf-8") as out:
        for i, pair in enumerate(one2one):
            sp1 = pair[0]
            sp2 = pair[1]
            score = expr_conservation_orthologs[i]
            out.write('\t'.join([sp1, sp2, str(score)])+'\n')


def write_ec_manyortho(results, out_many, out_skipped):

    """
    Writes scores to a tab-seprated file.

    Args:
        results (list): list with ec scores for all pairs of many-to-many / many-to-one orthologs
        out_many (numpy.array): name of output file to write scores to
        out_skipped (int): name of output file to write skipped multigenic families
    """

    nbmany, nbskipped = 0, 0
    with open(out_many, 'w', encoding='utf-8') as out,\
         open(out_skipped, 'w', encoding='utf-8') as out2:
        for i in results:
            if not i.startswith('skipped'):
                out.write(i)
                nbskipped += 1
            else:
                out2.write(i)
                nbmany += 1
    return nbmany, nbskipped


def icc(matrix_file_a, matrix_file_b, families_file, outprefix, max_combin=300,
        ncores=1, batch_size=10, ono2one_only=False, seed=None, bar_format=BAR_FORMAT):

    """
    Use the ICC approach to: (i) compute expression conservation scores for 1-to-1 orthologs, (ii)
    select best pairs (i.e. most conserved expression) for many-to-many / many-to-one / one-to-many
    orthologs.

    Args:
        matrix_file_a, matrix_file_b (str): Tab-delimited gene - cluster expression (FC)
                                            matrix files for species 1 and 2, respectively
        families_file (str): file with orthologous gene pairs (tab-delimited, 1 pair per line)
        outprefix (str): prefix for output files
        max_combin (int, optional): maximum accepted number or pairwise combinations,
                                    massively multigeneic families will be skipped.
        ncores (int, optional) number of cores to use
        batch_size (int, optional): size of batch for each parallel job
        ono2one_only (bool, optional): do not select best pairs for many-to-many / one-to-many,
                                       use only one-to-one orthologs
        seed (int, optional): random seed, use for reproducible results
        bar_format (str, optional): tqdm bar format, use None for tqdm default
    """

    #Load gene - cluster expression matrices (Fold-change)
    logger.info('Loading gene-cluster expression matrices')
    mat1 = load_cluster_expression_matrix(matrix_file_a)
    mat2 = load_cluster_expression_matrix(matrix_file_b)


    # Susbet matrices to retain 1-1 orthologs only
    logger.info('Loading gene orthologies')
    one2one, manyortho = load_orthologs(families_file, set(mat1.genes), set(mat2.genes))
    one2one_a = [mat1.genes[i[0]] for i in one2one]
    one2one_b = [mat2.genes[i[1]] for i in one2one]
    a = mat1.matrix[one2one_a,:]
    b = mat2.matrix[one2one_b,:]
    logger.info('Found %s one-to-one orthologs', len(one2one_a))

    if len(one2one_a) < 1000:
        logger.error('Too few one-to-one orthologs, please check your orthology file.')
        sys.exit(1)

    # Compute 1-1 orthologs co-expression conservation
    logger.info('Computing co-expression conservation for one-to-one orthologs')
    expr_conservation_orthologs = compute_expression_conservation(a, b)

    write_ec_one2one(one2one, expr_conservation_orthologs,
                     outprefix + '1-to-1-orthologs_correlation_scores.tsv')


    if not ono2one_only:

        #Select best pairs for many-many / one-many
        logger.info('Selecting best pair for many-to-many / one-to-many / many-to-one orthologs')

        try:

            pool = multiprocessing.Pool(ncores, init_worker) #, maxtasksperchild=200
            og_batches = [manyortho[i:i + batch_size] for i in range(0, len(manyortho), batch_size)]

            if seed:
                np.random.seed()

            idx_random = np.random.choice(a.shape[0], 1000, replace=False)

            a = a[idx_random, :]
            b = b[idx_random, :]

            mat1 = (mat1.matrix, mat1.genes)
            mat2 = (mat2.matrix, mat2.genes)

            jobs = [pool.apply_async(worker_add_gene_pairs, args=(batch, a, b, mat1, mat2,
                                                                  max_combin))
                                                                  for batch in og_batches]
            pool.close()

            if bar_format:
                bar_format = bar_format.replace('task', 'Homologs selection')

            prbar = tqdm.tqdm(jobs, colour='#595c79', bar_format=bar_format)
            prbar.unit = ""
            prbar.refresh()

            async_res = []
            for i, job in enumerate(prbar):
                async_res += job.get()
                if i == len(prbar) - 1:
                    prbar.unit = "[done]"
                    prbar.refresh()


        except KeyboardInterrupt:
            logger.info("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
            sys.exit(1)

    else:
        logger.info('Using 1-to-1 orthologs only')


    nmany, nskip = write_ec_manyortho(async_res, outprefix+'orthologs_many_correlation_scores.txt',
                                      outprefix+'skipped_ogs.txt')

    logger.info('Added %s many-to-many / one-to-many / many-to-one, %s multigenic families skipped.'
                ' Total retained genes for cross-species comparison: %s genes.', nmany, nskip,
                 nmany+nskip)
