"""
Module with functions to compute gene expression similarity between single-cell datasets of 2
species, using pearson weighted correlation and pre-computed gene weights.
"""

import sys

import logging
import pickle
import collections

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

import coloredlogs

from . import iterative_comparison_coexpression as icc

logger = logging.getLogger(__name__)
coloredlogs.install()

def load_ec_ortho_1to1(input_file):

    """
    Loads expression conservation scores of 1-1 orthologs.

    Args:
        input_file (str): path to tab-delimited input file (columns = gene sp1, gene sp2. score)

    Returns:
        tuple:
            two lists, 1) list of orthologous gene pairs and 2) list of corresponding
            co-expression conservation scores (weights).
    """
    scores, genes = [], []
    with open(input_file, 'r', encoding = "utf-8") as infile:
        for line in infile:
            gene1, gene2, score = line.strip().split('\t')
            genes.append(f'{gene1}+{gene2}')
            scores.append(float(score))
    return genes, scores


def load_ec_ortho_many(input_file, threshold=None):

    """
    Loads expression conservation scores of many-to-many / many-to-one orthologs.

    Args:
        input_file (str): path to input file (format = line 1 gene pairs, line 2 associated scores, 
                          line 3 gene pairs, line 4 scores, ... etc)
        threshold (float): if provided, instead of best pair retain all distinct pairs with
                           score > threshold

    Returns:
        tuple:
            two lists, 1) list of orthologous gene pairs and 2) list of corresponding
            co-expression conservation scores (weights).
    """
    res_genes = []
    res_ec = []
    # seen = set()
    with open(input_file, 'r', encoding = "utf-8") as infile:
        for i, line in enumerate(infile):
            if i%2 == 0:
                genes = line.strip().split()
            else:
                scores = line.strip().split()
                if threshold is None:
                    best_ind = np.argmax([float(i) for i in scores])
                    res_genes.append(genes[best_ind])
                    res_ec.append(float(scores[best_ind]))
                else:
                    ind_sorted = np.flip(np.argsort([float(i) for i in scores]))
                    for ind in ind_sorted:
                        if float(scores[ind]) <= threshold:
                            break

                        pair = genes[ind]
                        # g1, g2 = pair.split('+')
                        # if g1 not in seen and g2 not in seen: #TODO try even without this one
                        res_genes.append(pair)
                        res_ec.append(float(scores[ind]))
                            # seen.add(g1)
                            # seen.add(g2)
    return res_genes, res_ec


def filter_matrix(matrix, genedict, genes_to_keep):

    """
    Filters out expression matrix to retain only genes present in `genes_to_keep`.

    Args:
        matrix (numpy.array): expression matrix to filter
        genedict (dict): gene name (key) to matrix row indice (value)
        genes_to_keep (list): list of gene names to retain 

    Returns:
        numpy.array: filtered expression matrix
    """

    index_to_keep = [genedict[i] for i in genes_to_keep]
    res = matrix[index_to_keep,:]
    return res

def make_palette_for_broad(broad_file1, broad_file2, cell_types1, cell_types2):

    """
    Creates a color palette for broad annotations

    Args:
        broad_file1, broad_file2 (str): path to pickled dict giving cluster to broad correspondence
                                        for species1 and species2, respectively
        cell_types1, cell_types2 (list): list of clusters in species1 and species2, giving order in
                                         heatmap
    Returns:
        tuple:
            lists of broad annotations for species1 and species2 in the same order as clusters,
            and a dict with broad to color palette
    """

    #load dict that was pickled during input file loading
    with open(broad_file1, 'rb') as infile:
        broad1 = pickle.load(infile)

    with open(broad_file2, 'rb') as infile:
        broad2 = pickle.load(infile)

    broad1 = {key:broad1[key] for key in broad1 if key in cell_types1}
    broad2 = {key:broad2[key] for key in broad2 if key in cell_types2}
    broad_all = [broad1[c] for c in cell_types1] + [broad2[c] for c in cell_types2]

    unique_broad = sorted(list(set(broad_all)))

    colors = sns.color_palette("Set3", 12) + sns.color_palette("tab20", 20)
    palettedict = {}

    #hard-coding 'unknown' or 'other' to be gray and 'muscle' to be red
    if 'Unknown' in unique_broad or 'Other' in unique_broad:
        if 'Unknown' in unique_broad:
            palettedict['Unknown'] = colors[8]
            del colors[8]
        else:
            palettedict['Other'] = colors[8]
            del colors[8]

    if 'Muscle' in unique_broad:
        palettedict['Muscle'] = colors[3]
        del colors[3]

    if len(unique_broad) > len(colors):
        logger.warning('Found over 32 broad cell type classes (%s), broad annotation will not be '
                       'plotted.', len(unique_broad))
        return None

    for i, b in enumerate(unique_broad):
        if b not in palettedict:
            palettedict[b] = colors[i]
        else:
            i = i - 1

    return broad1, broad2, palettedict


def plot_and_save_out(result, cell_types1, cell_types2, outprefix, sp1='', sp2='',
                      threshold_for_plot=0, outformat='svg', broad_file1=None, broad_file2=None,
                      reorder='DiagKeep'):

    """
    Saves and plots heatmap showing expression comparisons between all clusters of sp1 and all
    clusters of sp2.

    Args:
        result (numpy.array): matrix with weighted pearson correlation for gene expression
                              in all pairs of clusters
        cell_types1, cell_types2 (list): list of clusters in rows (species2) and columns (species2)
        outprefix (str): prefix for output files
        sp1, sp2 (str, optional): (short)name of species 1 and species 2 (to label axes on matrix)
        threshold_for_plot (float, optional): set all result entry < threshold_for_plot to 0, to
                                              filter out possibly spurrious match, consider using
                                              pesci randomization to estimate threshold
        outformat (str, optional): format for saved figure (svg, png or pdf)
        broad_file1, broad_file2 (str, optional): path to pickle file with dict cluster to broad
                                                  annotations (file generated by pesci)
        reorder (str, optional): if not using broad annotation to order the heatmap, to alternatives
                                 are possible: 'Diagkeep' will put best 1-1 match on the Diagonal,
                                 'Clust' will use the order from hierarchical 
                                 clustering from seaborn.clustermap
    """
    df = pd.DataFrame(data=result[0:,0:], columns=[sp2+'|'+i for i in cell_types2],
                      index=[sp1+'|'+i for i in cell_types1])
    df.to_csv(f'{outprefix}correlation_scores_matrix.csv', sep='\t')

    #try to get cells to be ~square
    if len(cell_types2) / len(cell_types1) > 1.25:
        fsize = (8, 6)
        cratio = (0.006*3, 0.008*3)
    elif len(cell_types1) / len(cell_types2) > 1.25:
        fsize = (6, 8)
        cratio = (0.008*3, 0.006*3)
    else:
        fsize = (8, 8)
        cratio = (0.008*3, 0.008*3)

    result[result<threshold_for_plot] = 0

    row_colors, col_colors = None, None
    if broad_file1 and broad_file2:
        pal = make_palette_for_broad(broad_file1, broad_file2, cell_types1, cell_types2)
        if pal:
            broad1, broad2, paldict = pal
            c1 = [broad1[i] for i in cell_types1]
            c2 = [broad2[i] for i in cell_types2]

            sorter1 = sorted(range(len(broad1)), key=c1.__getitem__)
            cell_types1 = [list(cell_types1.keys())[i] for i in sorter1]
            result = result[sorter1,:]

            sorter2 = sorted(range(len(broad2)), key=c2.__getitem__)
            cell_types2 = [list(cell_types2.keys())[i] for i in sorter2]
            result = result[:,sorter2]

            row_colors = [paldict[broad1[i]] for i in cell_types1]
            col_colors = [paldict[broad2[i]] for i in cell_types2]

    elif reorder == 'DiagKeep':
        row_idx, col_idx = linear_sum_assignment(-result)
        all_row = range(len(result[:,0]))
        col_idx_final, row_idx_final = [], []
        for idx in all_row:
            if idx in row_idx:
                i = list(row_idx).index(idx)
                row_idx_final.append(idx)
                col_idx_final.append(col_idx[i])

        row_idx, col_idx = list(row_idx), list(col_idx)
        all_col = range(len(result[0,:]))

        row_idx_final = row_idx_final + [i for i in all_row if i not in row_idx_final]
        col_idx_final = col_idx_final + [i for i in all_col if i not in col_idx_final]
        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types2 = [list(cell_types2.keys())[i] for i in col_idx_final]
        cell_types1 = [list(cell_types1.keys())[i] for i in row_idx_final]


    elif reorder == 'Clust':
        #use clustermap to get a reordering
        clustergrid = sns.clustermap(result)

        row_idx_final = clustergrid.dendrogram_row.reordered_ind
        col_idx_final = clustergrid.dendrogram_col.reordered_ind

        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types2 = [list(cell_types2.keys())[i] for i in col_idx_final]
        cell_types1 = [list(cell_types1.keys())[i] for i in row_idx_final]
        plt.close('all')

    else:
        logger.error('%s is not a valid argument for the reorder parameter, '
                       'use DiagKeep or Clust', reorder)
        sys.exit(1)

    g = sns.clustermap(result, cmap='BuPu', annot=False, vmin=0, vmax=1, xticklabels=cell_types2,
                yticklabels=cell_types1, cbar_kws={'label': 'weighted\ncorrelation', "shrink": 0.1},
                row_cluster=False, col_cluster=False, dendrogram_ratio=0.01,
                figsize=fsize, cbar_pos=(1.05, 0.8, 0.02, 0.15), row_colors=row_colors,
                col_colors=col_colors, colors_ratio=cratio)

    if broad_file1 and broad_file2:
        handles = [Patch(facecolor=paldict[name]) for name in sorted(list(paldict.keys()))]
        plt.legend(handles, sorted(list(paldict.keys())), title='Broad annotation',
                   bbox_to_anchor=(1.22, 0.1), bbox_transform=plt.gcf().transFigure,
                   loc='lower right')

    ax = g.ax_heatmap
    ax.set_ylabel(sp1)
    ax.set_xlabel(sp2)
    plt.savefig(f'{outprefix}correlation_scores_matrix.{outformat}', bbox_inches='tight')
    plt.close('all')


def plot_expression_conservation(ec, n_ortho, outprefix):

    """
    Plots distribution of ec scores (histogram)

    Args:
        ec (numpy.array): vectoir with ec scores
        n_ortho (int): index at which many-to-many start in vector
        outprefix (str): prefix for output file
    """

    one2one_ec = [(i, 'one-to-one') for i in ec[:n_ortho]]
    manyortho_ec = [(i, 'other') for i in ec[n_ortho:]]
    data = manyortho_ec + one2one_ec
    df = pd.DataFrame.from_records(data, columns=['Expression conservation score',
                                                  'homologs type'])
    sns.histplot(data=df, x="Expression conservation score", hue="orthology type",
                 bins=41, multiple="stack")
    plt.tight_layout()
    plt.savefig(outprefix+'expression_correlation_scores.svg')
    plt.close('all')


def make_coexpressed_genes_table(result, mat1, mat2, ec, idx_ortho, outprefix, sp1='sp1', sp2='sp2',
                                 fc=1.5, wmin=0):

    """
    Searches for co-expressed marker genes for matching clusters in the pairwise cross species
    comparison.

    Args:
        result (numpy.array): matrix with weighted pearson correlation for gene expression
                              in all pairs of clusters
        mat1, mat2 (ExprMatrix2 namedtuple): expression matrix (numpy array), list of genes in rows
                                             and list of clusters in columns
        ec (numpy.array): expression conservation score vector (same order as in matrix rows)
        idx_ortho (int): index at which many-to-many orthologs start in matrix
        outprefix (str): prefix for output file
        sp1, sp2 (str, optional): (short)name for species 1 and species 2, respectively
        fc (float, optional): minimum fold change to be considered marker of a cluster
    """

    records = []
    for i, clus1 in enumerate(mat1.clusters):
        for j, clus2 in enumerate(mat2.clusters):
            if result[i, j] > wmin:
                genes_clus1_sp1 = np.where(mat1.matrix[:, i] > fc)[0]
                genes_clus2_sp2 = np.where(mat2.matrix[:, j] > fc)[0]
                ortho_idx = set(genes_clus1_sp1).intersection(set(genes_clus2_sp2))

                for k in ortho_idx:
                    g1, g2 = mat1.genes[k], mat2.genes[k]
                    fc1, fc2 = mat1.matrix[k, i], mat2.matrix[k, j]
                    ec_tmp = ec[k]
                    score_tmp = np.exp(np.log(fc1 + fc2)) * max(ec_tmp, 0.1)
                    o1 = 'y'
                    if k >= idx_ortho:
                        o1 = 'n'
                    records.append([sp1+'|'+clus1, sp2+'|'+clus2, g1, g2, mat1.matrix[k, i],
                                    mat2.matrix[k, j], ec[k], score_tmp, o1])

    df = pd.DataFrame.from_records(records, columns=[f'{sp1}_cell_cluster', f'{sp2}_cell_cluster',
                                                     f'{sp1}_gene', f'{sp2}_gene',
                                                     f'{sp1}_Expression (Fold-change)', 
                                                     f'{sp2}_Expression (Fold-change)', 'EC_score',
                                                     'tmp_score', 'One-to-one-ortholog (yes/no)'])
    df.sort_values([f'{sp1}_cell_cluster', f'{sp2}_cell_cluster', 'tmp_score'], ascending=False,
                    inplace=True)
    df.drop(columns=['tmp_score'], inplace=True)
    df.to_csv(f'{outprefix}gene_coexpression_table.csv', sep='\t', index=False)


def compare(matrix_a, matrix_b, outprefix, sp1='sp1', sp2='sp2', random_id='',
            threshold_for_plot=0, outformat='svg', min_fc=1.5, broad_file1=None,
            broad_file2=None, many_threshold=None):

    """
    Compares gene expression across all pairs clusters species 1 - clusters species 2, using
    weighted pearson correlation. Resulting comparison heatmap is saved as figure and tab-delimited
    text file. Also reports co-expressed gene markers for matching clusters.

    Args:
        matrix_a, matrix_b (str): path to input file with gene - cluster normalized expression
                                  matrices (Fold-change)
        outprefix (str): prefix for output files
        sp1, sp2 (str, optional): (short)name for species 1 and species 2, respectively
        random_id (str, optional): set param to specify that orthologs have been randomized
                                   (results will be saved in a random_id subfolder,
                                   so consider setting random_id=random1 or random2,... etc)
        threshold_for_plot (str, optional): threshold to plot matches in comparison matrix & search
                                            for co-expressed marker gene pairs
        outformat (str, optional): format for saved figure (svg, png or pdf)
        min_fc (float, optional): minimum fold change to be considered marker of a cluster
        broad_file1, broad_file2 (str, optional): path to pickle file with dict cluster to broad
                                                  annotations (file generated by pesci)
        many_threshold (float): if provided, instead of best pair retain all distinct pairs with
                                score > threshold
    """

    #load expression matrices
    mat1 = icc.load_cluster_expression_matrix(matrix_a)
    mat2 = icc.load_cluster_expression_matrix(matrix_b)

    #load weight to compute correlations between cell types
    outprefix1 = outprefix +'files/'
    if random_id:
        outprefix1 = outprefix + random_id

    ortho, orthologs_ec = load_ec_ortho_1to1(outprefix1+sp1+'-'+sp2+\
                                             '_1-to-1-orthologs_correlation_scores.tsv')
    para, paralogs_ec = load_ec_ortho_many(outprefix1+sp1+'-'+sp2+\
                                           '_orthologs_many_correlation_scores.txt',
                                           threshold=many_threshold)

    ec = np.array(orthologs_ec + paralogs_ec)
    assert len(ec) == (len(para) + len(ortho))

    #filter matrix to retain 1-1 orthologs and selected best for many-to-many / one-to-many
    genes1_ok = [i.split('+')[0] for i in ortho+para]
    genes2_ok = [i.split('+')[1] for i in ortho+para]

    mat1_ok = filter_matrix(mat1.matrix, mat1.genes, genes1_ok)
    mat2_ok = filter_matrix(mat2.matrix, mat2.genes, genes2_ok)

    ExprMatrix2 = collections.namedtuple('ExprMatrix2', ['matrix', 'genes', 'clusters'])
    mat1 = ExprMatrix2(mat1_ok, genes1_ok, mat1.cells)
    mat2 = ExprMatrix2(mat2_ok, genes2_ok, mat2.cells)

    #TODO: remove this at the end when we are sure it's working as expected (unit tests needed)
    # df = pd.DataFrame(data=mat1_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para],
    #                   columns=[sp1+'|'+i for i in cell_types1])
    # df.to_csv(outprefix + 'file/' + Path(matrix_a).stem + '_orthologs.tsv', sep='\t')

    # df = pd.DataFrame(data=mat2_ok[0:,0:], index=[i.split('+')[0] for i in ortho+para],
    #                   columns=[sp2+'|'+i for i in cell_types2])
    # df.to_csv(outprefix + 'file/' + Path(matrix_b).stem + '_orthologs.tsv', sep='\t')

    is_one2one = np.concatenate([np.ones(len(ortho)), np.zeros(len(para))])
    df = pd.DataFrame([genes1_ok, genes2_ok, ec, is_one2one]).T
    df.columns = ['genes sp1', 'genes sp2', 'expression conservation', 'is_one2one']
    df.to_csv(outprefix+random_id+sp1+'-'+sp2+'_'+'expression_conservation_scores.csv',
              sep='\t', index=False)
    #compute weighted correlations between cell types of sp1 and cell types of sp2
    result = icc.column_wise_wcorr_einsum(mat1_ok, mat2_ok, ec)

    #get best match to print info
    maxind = np.unravel_index(result.argmax(), result.shape)
    maxval = round(result[maxind], 3)
    maxc1 = list(mat1.clusters.keys())[int(maxind[0])]
    maxc2 = list(mat2.clusters.keys())[int(maxind[1])]

    with open(outprefix1+f'{sp1}-{sp2}_best_match.txt', 'w', encoding='utf-8') as out:
        out.write('Max expression correlation value\tCluster sp1\tCluster sp2\n')
        out.write('\t'.join([str(maxval), maxc1, maxc2])+'\n')

    logger.info('Max expression correlation value: %s (%s|%s - %s|%s)', maxval, sp1, maxc1, sp2,
                                                                        maxc2)

    logger.info('Saving outputs and plotting correlation matrix')

    plot_and_save_out(result, mat1.clusters, mat2.clusters,
                               outprefix+random_id+sp1+'-'+sp2+'_', sp1, sp2,
                               threshold_for_plot=threshold_for_plot, outformat=outformat,
                               broad_file1=broad_file1, broad_file2=broad_file2)

    logger.info('Searching for co-expressed gene pairs')

    make_coexpressed_genes_table(result, mat1, mat2, ec, len(ortho),
                                 outprefix+random_id+sp1+'-'+sp2+'_', sp1, sp2, fc=min_fc,
                                 wmin=threshold_for_plot)


    # out_genes = outprefix+f'ortho_and_para.txt'
    # homologs = ortho+para
    # with open(out_genes, 'w') as out:
    #     for pair in homologs:
    #         out.write(pair+'\n')
