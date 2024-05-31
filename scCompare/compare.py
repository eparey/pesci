import sys
import logging

import multiprocessing
import traceback
import signal

import tqdm

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment

import seaborn as sns
import matplotlib.pyplot as plt

import coloredlogs

from . import iterative_comparison_coexpression as icc

logger = logging.getLogger(__name__)
coloredlogs.install()

def load_ec_ortho(input_file):
    scores, genes = [], []
    with open(input_file, 'r') as infile:
        for line in infile:
                gene1, gene2, score = line.strip().split('\t')
                genes.append(f'{gene1}+{gene2}')
                scores.append(float(score))
    return genes, scores


def parse_ec_para(input_file):
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
    return res_genes, res_ec


def filter_matrix(matrix, genedict, genes_to_keep):
    index_to_keep = [genedict[i] for i in genes_to_keep]
    res = matrix[index_to_keep,:]
    return res


def plot_and_save_out(result, cell_types2, cell_types1, outprefix, sp1='', sp2='', suffix='',
                      reorder='None'):

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

    df = pd.DataFrame(data=result[0:,0:], index=[sp2+'|'+i for i in cell_types2],
                      columns=[sp1+'|'+i for i in cell_types1])
    df.to_csv(outprefix+f'correlation_scores_matrix{suffix}.csv', sep='\t')

    plt.figure(figsize=(6, 6))
    sns.heatmap(result, cmap='BuPu', annot=False, vmin=0, vmax=1, xticklabels=cell_types1,
                yticklabels=cell_types2, cbar_kws={'label': 'weighted correlation', "shrink": 0.5},
                square=True)
    plt.ylabel(sp1)
    plt.xlabel(sp2)
    # plt.tight_layout()
    plt.savefig(outprefix+f'correlation_scores_matrix{suffix}.svg', bbox_inches='tight')
    plt.close('all')


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def plot_expression_conservation(ec, n_ortho, outprefix, suffix=''):

    orthologs_ec = [(i, 'orthologs') for i in ec[:n_ortho]]
    paralogs_ec = [(i, 'paralogs') for i in ec[n_ortho:]]
    data = paralogs_ec + orthologs_ec
    df = pd.DataFrame.from_records(data, columns=['Expression conservation score', 'homologs type'])
    sns.histplot(data=df, x="Expression conservation score", hue="homologs type",
                 bins=41, multiple="stack")
    plt.tight_layout()
    plt.savefig(outprefix+f'expression_correlation_scores{suffix}.svg')
    plt.close('all')


def compare_main(matrix_a, matrix_b, outprefix, ncores, sp1='', sp2=''):
    #load expression matrices
    mat1, genes1, cell_types1 = icc.parse_matrix(matrix_a)
    mat2, genes2, cell_types2 = icc.parse_matrix(matrix_b)

    #load weight to compute correlations between cell types
    ortho_ec = outprefix + 'one-to-one-orthologs_correlation_scores.csv'
    ortho, orthologs_ec = load_ec_ortho(ortho_ec)
    para_ec = outprefix + 'paralogs_correlation_scores.txt'
    para, paralogs_ec = parse_ec_para(para_ec)

    ec = np.array(orthologs_ec + paralogs_ec)
    assert len(ec) == (len(para) + len(ortho))
    # logger.info(f'Loaded gene weigths ({len(ec)} genes)')

    #filter matrix to retain 1-1 orthologs and selected best paralogs
    mat1_ok = filter_matrix(mat1, genes1, [i.split('+')[0] for i in ortho+para])
    mat2_ok = filter_matrix(mat2, genes2, [i.split('+')[1] for i in ortho+para])

    df = pd.DataFrame(data=mat1_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para],
                      columns=[sp1+'|'+i for i in cell_types1])
    df.to_csv(outprefix+f'matrix_{sp1}.csv', sep='\t')

    df = pd.DataFrame(data=mat2_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para],
                      columns=[sp2+'|'+i for i in cell_types2])
    df.to_csv(outprefix+f'matrix_{sp2}.csv', sep='\t')

    df = pd.DataFrame(data=ec, index=[i.split('+')[0] for i in ortho+para], columns=['ec'])
    df.to_csv(outprefix+'ec_vector.csv', sep='\t')

    #compute weighted correlations between cell types of sp1 and cell types of sp2
    result = icc.column_wise_wcorr_einsum(mat1_ok, mat2_ok, ec)

    logger.info('Saving outputs and plotting correlation matrix')

    plot_and_save_out(result, cell_types1, cell_types2, outprefix, sp1, sp2)

    plot_expression_conservation(ec, len(ortho), outprefix)

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
