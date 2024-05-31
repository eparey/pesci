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

logger = logging.getLogger(__name__)
coloredlogs.install()


def wcov(v1, v2, w):
    return np.sum(w * (v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w))) / np.sum(w)
    # return np.average((v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w)), weights=w)

def wcorr(v1, v2, w, get_markers=False):
    res = wcov(v1, v2, w) / np.sqrt(wcov(v1, v1, w) * wcov(v2, v2, w))
    if get_markers:
        all_contrib = w * (v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w))
        n_tops = 10
        idx_best = np.argpartition(all_contrib, -n_tops)[-n_tops:]
        v1_best = v1[idx_best]
        v2_best = v2[idx_best]
        w_best = w[idx_best]
        return res, [idx_best, v1_best, v2_best, w_best, all_contrib[idx_best]]

    return res


def init_ec(mat1, mat2, n):
    ec = []
    for i in range(n):
        ec.append(max(np.corrcoef(mat1[i,:], mat2[i,:])[0,1], 0))
    return np.array(ec)


def init_ec_optimized_einsum(mat1, mat2, n):

    mat1_m = mat1 - np.einsum("ij->i", mat1, optimize='optimal') / np.double(n)
    mat2_m = mat2 - np.einsum("ij->i", mat2, optimize='optimal') / np.double(n)

    cov = np.einsum("ij,ij->j", mat1_m, mat2_m, optimize='optimal')

    var1 = np.einsum("ij,ij->j", mat1_m, mat1_m, optimize='optimal')
    var2 = np.einsum("ij,ij->j", mat2_m, mat2_m, optimize='optimal')

    varprod = np.einsum("i,i->i", var1, var2, optimize='optimal')

    ec = cov / np.sqrt(varprod)
    return ec.clip(min=0)


def update_ec(mat1, mat2, prev_ec, n):
    ec = []
    for i in range(n):
        ec.append(max(wcorr(mat1[i,:], mat2[i,:], prev_ec), 0))
    return np.array(ec)


def update_ec_optimized_einsum(mat1, mat2, prev_ec, n):

    wmat = np.tile(prev_ec, (n, 1))
    mat1_m = mat1 - np.einsum("ji,ij->j", wmat, mat1, optimize='optimal') / sum(prev_ec)
    mat2_m = mat2 - np.einsum("ji,ij->j", wmat, mat2, optimize='optimal') / sum(prev_ec)

    cov = np.einsum("ij,ij->ij", mat1_m, mat2_m, optimize='optimal')
    wcov12 = np.einsum("ji,ij->j", wmat, cov, optimize='optimal') / sum(prev_ec)

    cov1 = np.einsum("ij,ij->ij", mat1_m, mat1_m, optimize='optimal')
    wcov1 = np.einsum("ji,ij->j", wmat, cov1, optimize='optimal') / sum(prev_ec)

    cov2 = np.einsum("ij,ij->ij", mat2_m, mat2_m, optimize='optimal')
    wcov2 = np.einsum("ji,ij->j", wmat, cov2, optimize='optimal') / sum(prev_ec)

    prod = np.einsum("i,i->i", wcov1, wcov2, optimize='optimal')
    res = wcov12 / np.sqrt(prod)
    res[np.isnan(res)] = 0
    return res.clip(min=0)


def update_ec_until_convergence(mat1, mat2, ec, n, max_iter=100):
    i = 0
    diff = np.Inf
    while i < max_iter and diff > 0.05:
        ec_new = update_ec_optimized_einsum(mat1, mat2, ec, n)
        diff = np.sum((ec - ec_new)**2)
        ec = ec_new
        i += 1

    if i == max_iter:
        logger.warning('Warning, ICC convergence not reached, consider increasing max_iter.')

    return ec


def compute_expression_conservation(matrix_sp_a, matrix_sp_b):

    # Transform each expression matrix into a correlation co-expression matrix
    co_a = np.corrcoef(matrix_sp_a)
    co_b = np.corrcoef(matrix_sp_b)

    n = len(co_a[:,0])

    #For each row in a, compute correlation with the corresponding row in b to init EC
    expr_conservation = init_ec_optimized_einsum(co_a, co_b, n)

    #Update EC using weithed correlation until convergence or max iter is reached
    expr_conservation_final = update_ec_until_convergence(co_a, co_b, expr_conservation, n)

    return expr_conservation_final


def load_orthologs(input_file, genes_sp_a, genes_sp_b):
    edges = []
    with open(input_file, 'r') as infile:
        for line in infile:
            g1, g2 = line.strip().split('\t')
            if g1 in genes_sp_a and g2 in genes_sp_b:
                edges.append((g1, g2))
            elif g2 in genes_sp_a and g1 in genes_sp_b:
                edges.append((g2, g1))
    # print(f"{len(edges)} edges")
    homologs_graph = nx.Graph()
    homologs_graph.add_edges_from(edges)


    components = nx.connected_components(homologs_graph)
    # print('Connected components extracted')
    orthologs, paralogs = [], []
    for component in components:
        genes = set(component)
        if len(component) == 2:
            g1, g2 = genes
            if g1 in genes_sp_b:
                g1, g2 = g2, g1
            orthologs.append((g1, g2))
        else:
            tmp_genes_a = genes.intersection(genes_sp_a)
            tmp_genes_b = genes.intersection(genes_sp_b)
            paralogs.append((tuple(tmp_genes_a), tuple(tmp_genes_b)))
    return orthologs, paralogs

def parse_matrix(input_file):

    #load data
    df = pd.read_csv(input_file, sep="\t", header=0, index_col=0)

    #remove constant rows
    unique_values = df.nunique(axis=1)
    constant_rows = unique_values[unique_values == 1].index.tolist()
    df = df.drop(constant_rows)

    #get gene names
    genes = df.index.values.tolist()
    genes = {genes[i]:i for i in range(len(genes))}

    #get cell types
    cell_types = df.columns.tolist()

    return np.array(df), genes, cell_types


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    #create a unique seed generator for each pool
    #TDOD: check this works as expected
    global rng
    rng = np.random.default_rng()


def worker_best_paralog(paralogs, a, b, mat1, mat2, mat1_genes, mat2_genes, max_combin):

    try:

        res = []
        for og in paralogs:

            #list all possible pairings (create list of tuples)
            combin = [(i,j) for (i,j) in list(itertools.product(og[0], og[1])) if i in mat1_genes and j in mat2_genes]

            if len(combin) > max_combin:
                res.append(f'skipped {og} - {len(og[0])} sp1 genes - {len(og[1])} sp2 genes - {len(combin)} combinations\n')

            else:
                #add all 'a' paralogs to matrix a (with duplicates to accomodate all possible pairings)
                paralogs_a = [mat1_genes[i[0]] for i in combin]

                idx_random = np.random.choice(a.shape[0], 1000, replace=False)

                small_a = a[idx_random, :] #subset a to reduce time (1000 random orthologs), or we only do once for all selected paralogs??
                a_w_paralogs_tmp = np.concatenate((small_a, mat1[paralogs_a,:]), axis=0)

                #add all 'b' paralogs to matrix b (with duplicates to accomodate all possible pairings)
                paralogs_b = [mat2_genes[i[1]] for i in combin]

                small_b = b[idx_random, :] #subset b to reduce time (1000 random orthologs)
                b_w_paralogs_tmp = np.concatenate((small_b, mat2[paralogs_b,:]), axis=0)

                #compute expresseion conservation with all possible 1-1 groupings of paralogs
                expr_conservation_w_paralogs_current = compute_expression_conservation(a_w_paralogs_tmp, b_w_paralogs_tmp)

                tmp = '\t'.join(['+'.join(comb) for comb in combin])+ '\n' + '\t'.join([str(i) for i in expr_conservation_w_paralogs_current[1000:]])+'\n'
                res.append(tmp)
        return res

    except Exception:
        traceback.print_exc()
        raise



def select_paralogs_main(matrix_file_a, matrix_file_b, families_file, outprefix, max_combin=2500,
                         ncores=1, batch_size=5, noparalogs=False):

    logger.info('Loading gene-cluster expression matrices')

    mat1, mat1_genes, _ = parse_matrix(matrix_file_a)

    mat2, mat2_genes, _ = parse_matrix(matrix_file_b)


    # Susbet the matrix to retain only 1-1 orthologs
    logger.info('Parsing gene families')
    orthologs, paralogs = load_orthologs(families_file, set(mat1_genes), set(mat2_genes))

    orthologs_a = [mat1_genes[i[0]] for i in orthologs]
    orthologs_b = [mat2_genes[i[1]] for i in orthologs]
    a = mat1[orthologs_a,:]
    b = mat2[orthologs_b,:]
    n = len(orthologs_b)
    logger.info(f'Found {n} one-to-one orthologs')

    if n < 1000:
        logger.error('Too few one-to-one orthologs, please check your orthology file.')
        sys.exit(1)

    # Compute orthologs co-expression conservation
    logger.info('Computing expression conservation of one-to-one orthologs')
    expr_conservation_orthologs = compute_expression_conservation(a, b)

    out_ortho = outprefix + 'orthologs_correlation_scores.txt'
    with open(out_ortho, 'w') as out:
        names_a = [i[0] for i in orthologs]
        names_b = [i[1] for i in orthologs]
        out.write('\t'.join([names_a[i]+'+'+names_b[i] for i in range(n)])+'\n')
        out.write('\t'.join([str(i) for i in expr_conservation_orthologs])+'\n')

    out_ortho = outprefix + 'orthologs_correlation_scores.csv'
    with open(out_ortho, 'w') as out:
        for i, pair in enumerate(orthologs):
            sp1 = pair[0]
            sp2 = pair[1]
            score = expr_conservation_orthologs[i]
            out.write('\t'.join([sp1, sp2, str(score)])+'\n')

    logger.info('Selecting best homologs for remaining gene families')
    async_res = []
    if not noparalogs:
        try:

            pool = multiprocessing.Pool(ncores, init_worker) #, maxtasksperchild=200
            og_batches = [paralogs[i:i + batch_size] for i in range(0, len(paralogs), batch_size)]
            jobs = [pool.apply_async(worker_best_paralog, args=(batch, a, b, mat1, mat2, mat1_genes, mat2_genes, max_combin)) for batch in og_batches]
            pool.close()
            # pool.join()
            pbar = tqdm.tqdm(jobs, colour='#595c79', bar_format='{percentage:3.0f}% |{bar:50}| Homologs selection \x1B[1;32m{unit}')
            pbar.unit = ""
            pbar.refresh()
            for i, job in enumerate(pbar):
                async_res += job.get()
                # print(len(async_res))
                if i == len(pbar) - 1:
                    pbar.unit = "[done]"
                    pbar.refresh()


        except KeyboardInterrupt:
            logger.info("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
            sys.exit(1)

    out_para = outprefix + 'paralogs_correlation_scores.txt'
    out_skipped = outprefix + 'skipped_ogs.txt'
    h = 0
    k = 0
    with open(out_para, 'w') as out, open(out_skipped, 'w') as out2:
        for i in async_res:
            if not i.startswith('skipped'):
                out.write(i)
                h += 1
            else:
                out2.write(i)
                k += 1

    logger.info(f'Added {h} homologs, {k} multigenic families skipped. Total retained genes for cross-species comparison: {h+n} genes.')
