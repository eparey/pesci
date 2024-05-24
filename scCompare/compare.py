import sys

import collections
import itertools
import argparse
import tqdm

import multiprocessing
import traceback
import signal

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment

import seaborn as sns
import matplotlib.pyplot as plt

# import pingouin as pg

# import qnorm

from . import select_paralogs as sp


def parse_ec_ortho(input_file, format='single'):
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            if i%2 == 0:
                genes = line.strip().split()
            else:
                scores = line.strip().split()
    return genes, [float(i) for i in scores]


def parse_ec_para(input_file, format='single'):
    res_genes = []
    res_ec = []
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            if i%2 == 0:
                genes = line.strip().split()
            else:
                scores = line.strip().split()
                best_ind = np.argmax([float(i) for i in scores])
                res_genes.append(genes[best_ind])
                res_ec.append(float(scores[best_ind]))
                # if len(scores) > 3:
                #     scores.remove(scores[best_ind])
                #     best_ind = np.argmax([float(i) for i in scores])
                #     res_genes.append(genes[best_ind])
                #     res_ec.append(float(scores[best_ind]))
                # scores = line.strip().split()
                # for i in range(len(scores)):
                #     res_genes.append(genes[i])
                #     res_ec.append(float(scores[i]))
    return res_genes, res_ec


def filter_matrix(matrix, genedict, genes_to_keep):
    index_to_keep = [genedict[i] for i in genes_to_keep]
    res = matrix[index_to_keep,:]
    return res


def plot_and_save_out(result, cell_types1, cell_types2, outprefix, sp1='', sp2='', suffix='', reorder='None'):

    # #This does not work
    if reorder == 'Diag':
        row_idx_final, col_idx_final = linear_sum_assignment(-result) 
        all_row = list(range(len(result[:,0])))
        col_idx_final, row_idx_final = list(col_idx_final), list(row_idx_final)
        # for idx in all_row:
        #     if idx in row_idx:
        #         i = list(row_idx).index(idx)
        #         row_idx_final.append(idx)
        #         col_idx_final.append(col_idx[i])

        # print(row_idx, col_idx)
        # row_idx, col_idx = list(row_idx), list(col_idx)
        all_col = list(range(len(result[0,:])))
        
        row_idx_final = row_idx_final + [i for i in all_row if i not in row_idx_final]
        col_idx_final = col_idx_final + [i for i in all_col if i not in col_idx_final]
        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types1 = [cell_types1[i] for i in col_idx_final]
        cell_types2 = [cell_types2[i] for i in row_idx_final]

    if reorder == 'DiagKeep':
        row_idx, col_idx = linear_sum_assignment(-result) 
        all_row = range(len(result[:,0]))
        col_idx_final, row_idx_final = [], []
        for idx in all_row:
            if idx in row_idx:
                i = list(row_idx).index(idx)
                row_idx_final.append(idx)
                col_idx_final.append(col_idx[i])

        print(row_idx, col_idx)
        row_idx, col_idx = list(row_idx), list(col_idx)
        all_col = range(len(result[0,:]))
        
        row_idx_final = row_idx_final + [i for i in all_row if i not in row_idx_final]
        col_idx_final = col_idx_final + [i for i in all_col if i not in col_idx_final]
        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types1 = [cell_types1[i] for i in col_idx_final]
        cell_types2 = [cell_types2[i] for i in row_idx_final]

    if reorder == 'Clust':
        #use clustermap to get a reordering
        clustergrid = sns.clustermap(result)
    
        row_idx_final = clustergrid.dendrogram_row.reordered_ind
        col_idx_final = clustergrid.dendrogram_col.reordered_ind

        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types1 = [cell_types1[i] for i in col_idx_final]
        cell_types2 = [cell_types2[i] for i in row_idx_final]
        plt.close('all')

    df = pd.DataFrame(data=result[0:,0:], index=[sp2+'|'+i for i in cell_types2], columns=[sp1+'|'+i for i in cell_types1])
    df.to_csv(outprefix+f'correlation_scores_matrix{suffix}.csv', sep='\t')
    # np.savetxt(outprefix+f'correlation_scores_matrix{suffix}.txt', result)

    plt.figure(figsize=(10, 10))
    sns.heatmap(result, cmap='BuPu', annot=False, vmin=0, vmax=1, xticklabels=cell_types1, yticklabels=cell_types2, cbar_kws={'label': 'weighted correlation'})
    plt.ylabel(sp2)
    plt.xlabel(sp1)
    # plt.tight_layout()
    plt.savefig(outprefix+f'correlation_scores_matrix{suffix}.svg', bbox_inches='tight')
    plt.close('all')


def weighted_correlation_celltypes_pairs(cell_types1, cell_types2, mat1, mat2, ec, reoptimize_ec=False, mark=False):
    #compute weighted correlations between cell types of sp1 and cell types of sp2
    result = np.zeros((len(cell_types2), len(cell_types1)))

    # qnormed_mats = qnorm.quantile_normalize(np.concatenate((mat1, mat2), axis=1))
    # print(qnormed_mats)
    # print(qnormed_mats.shape)

    # idx_mat1 = [i for i in range(len(cell_types1))]

    # idx_mat2 = [i+len(cell_types1) for i in range(len(cell_types2))]

    # mat1 = qnormed_mats[:,idx_mat1]
    # mat2 = qnormed_mats[:,idx_mat2]
    # print(mat1.shape)

    if reoptimize_ec:
        ec = sp.compute_expression_conservation(mat1, mat2)

    markers = {}
    for i in range(len(cell_types2)):
        for j in range(len(cell_types1)):
            
            if mark:
                wcorr, top500 = sp.wcorr(mat2[:,i], mat1[:,j], ec, get_markers=mark)
                markers[(i, j)] = top500
            else:
                wcorr = sp.wcorr(mat2[:,i], mat1[:,j], ec)
            result[i,j] = max(wcorr, 0) #Loop here to remove in optimization
    if mark:
        return result, markers, ec
    else:
        return result, ec


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def worker_random_corr(mat1_ok, mat2_ok, cell_types1, cell_types2, i, outprefix, sp1, sp2):

    try:
        tmp = mat2_ok
        np.random.shuffle(tmp)
        expr_conservation = sp.compute_expression_conservation(mat1_ok, tmp)
        rand_corr = weighted_correlation_celltypes_pairs(cell_types1, cell_types2, mat1_ok, tmp, expr_conservation)
        np.savetxt(outprefix+f'random/correlation_scores_matrix_{i}.txt', rand_corr)
        return rand_corr

    except Exception:
        traceback.print_exc()
        raise


# def significant_matches(result, cell_types1, cell_types2, outprefix, all_random_scores, sp1='', sp2=''):
#     tests = []
#     pvals = []
#     max_score = 0
#     for i in range(len(cell_types2)):
#         for j in range(len(cell_types1)):
#             max_score_current = max(all_random_scores)
#             if max_score_current > max_score:
#                 max_score = max_score_current

#             score = result[i, j]
#             test = (cell_types1[j]+'_'+sp1, cell_types2[i]+'_'+sp2)
#             pval = len([s for s in all_random_scores[:,j] if s > score]) / len(all_random_scores)
#             # if pval == 0:
#             #     pval = 1 / len(all_random_scores)
#             tests.append(test)
#             pvals.append(pval)
#     rej, padj = pg.multicomp(pvals, alpha=0.05, method='fdr_by')
#     for i, p in enumerate(padj):
#         if rej[i]:
#             print(tests[i], p)
#     print(max_score)

def plot_expression_conservation(ec, n_ortho, outprefix, suffix=''):

    orthologs_ec = [(i, 'orthologs') for i in ec[:n_ortho]]
    paralogs_ec = [(i, 'paralogs') for i in ec[n_ortho:]]
    data = paralogs_ec + orthologs_ec
    df = pd.DataFrame.from_records(data, columns=['Expression conservation score', 'homologs type'])
    sns.histplot(data=df, x="Expression conservation score", hue="homologs type", bins=41, multiple="stack")
    plt.tight_layout()
    plt.savefig(outprefix+f'expression_correlation_scores{suffix}.svg')
    plt.close('all')


def compare_main(matrix_a, matrix_b, outprefix, ncores, sp1='', sp2='', n=100, reoptimize_ec=False):
    #load expression matrices
    mat1, genes1, cell_types1 = sp.parse_matrix(matrix_a)
    mat2, genes2, cell_types2 = sp.parse_matrix(matrix_b)

    #load weight to compute correlations between cell types
    ortho_ec = outprefix + 'orthologs_correlation_scores.txt'
    ortho, orthologs_ec = parse_ec_ortho(ortho_ec, format='single')
    para_ec = outprefix + 'paralogs_correlation_scores.txt'
    para, paralogs_ec = parse_ec_para(para_ec, format='single')

    print(len(ortho))
    # print(len({i.split('+')[0] for i in ortho}))
    # print(len({i.split('+')[1] for i in ortho}))

    print(len(para))
    # print(len({i.split('+')[0] for i in para}))
    # print(len({i.split('+')[1] for i in para}))

    # print([item for item, count in collections.Counter([i.split('+')[0] for i in para]).items() if count > 1])
    # print([item for item, count in collections.Counter([i.split('+')[1] for i in para]).items() if count > 1])

    # print(len(ortho+para))
    # print(len({i.split('+')[0] for i in ortho+para}))
    # print(len({i.split('+')[1] for i in ortho+para}))

    ec = np.array(orthologs_ec + paralogs_ec)
    print(len(ec))

    #filter matrix to retain 1-1 orthologs and selected best paralogs
    mat1_ok = filter_matrix(mat1, genes1, [i.split('+')[0] for i in ortho+para])
    mat2_ok = filter_matrix(mat2, genes2, [i.split('+')[1] for i in ortho+para])

    df = pd.DataFrame(data=mat1_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para], columns=[sp1+'|'+i for i in cell_types1])
    df.to_csv(outprefix+f'matrix_{sp1}.csv', sep='\t')

    df = pd.DataFrame(data=mat2_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para], columns=[sp2+'|'+i for i in cell_types2])
    df.to_csv(outprefix+f'matrix_{sp2}.csv', sep='\t')

    df = pd.DataFrame(data=ec, index=[i.split('+')[0] for i in ortho+para], columns=['ec'])
    df.to_csv(outprefix+f'ec_vector.csv', sep='\t')

    #compute weighted correlations between cell types of sp1 and cell types of sp2
    result, tops, ec_final = weighted_correlation_celltypes_pairs(cell_types1, cell_types2, mat1_ok, mat2_ok, ec, mark=True, reoptimize_ec=reoptimize_ec)

    plot_and_save_out(result, cell_types1, cell_types2, outprefix, sp1, sp2)

    plot_expression_conservation(ec_final, len(ortho), outprefix)

    # out_genes = outprefix+f'ortho_and_para.txt'
    # homologs = ortho+para
    # with open(out_genes, 'w') as out:
    #     for pair in homologs:
    #         out.write(pair+'\n')

    # out_genes = outprefix+f'top_genes.txt'
    # with open(out_genes, 'w') as out:
    #     for i in range(len(cell_types2)):
    #         for j in range(len(cell_types1)):
    #             top = tops[i, j]
    #             out.write(f'{cell_types1[j]} {sp1} - {cell_types2[i]} {sp2}\n')
    #             for k in range(len(top[0])):
    #                 idx, v1 , v2, w, score = top[0][k], top[1][k], top[2][k], top[3][k], top[4][k]
    #                 homologs = ortho+para
    #                 # print(homologs)
    #                 genes = homologs[idx]
    #                 out.write(f'{genes} {v1} {v2} {w} {score}\n')
    #             out.write(f'//\n')


    # #100 random correlation matrices
    # try:
    #     res_random = []
    #     pool = multiprocessing.Pool(ncores, init_worker) #, maxtasksperchild=200
    #     jobs = [pool.apply_async(worker_random_corr, args=(mat1_ok, mat2_ok, cell_types1, cell_types2, i, outprefix, sp1, sp2)) for i in range(n)]
    #     pool.close()
    #     # pool.join()
        
    #     for job in tqdm.tqdm(jobs, colour='#FF79C6'):
    #         res_random.append(job.get())


    # except KeyboardInterrupt:
    #     print("Caught KeyboardInterrupt, terminating workers")
    #     pool.terminate()
    #     pool.join()
    #     sys.exit(1)

    # significant_matches(result, cell_types1, cell_types2, outprefix, np.concatenate(res_random, axis=1), sp1, sp2)
