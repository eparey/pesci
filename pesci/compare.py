"""
Module with functions to compute gene expression similarity between single-cell datasets of 2
species, using pearson weighted correlation and pre-computed gene weights.
"""

import logging
import pickle
import collections

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment

from skimage.filters import threshold_otsu

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe
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

    #FIXME allow to provide user-defined colors
    if 'Meiotic' not in unique_broad and 'Peptidergic' in unique_broad:
        del colors[7]

    if 'Muscle' in unique_broad:
        palettedict['Muscle'] = colors[3]
        del colors[3]

    if len(unique_broad) > len(colors):
        logger.warning('Found over 32 broad cell type classes (%s), broad annotations will not be '
                       'plotted.', len(unique_broad))
        return None

    for i, b in enumerate(unique_broad):
        if b not in palettedict:
            palettedict[b] = colors[i]
        else:
            i = i - 1

    return broad1, broad2, palettedict


def plot_and_save_out(result, cell_types1, cell_types2, outprefix, sp1='', sp2='',
                      outformat='svg', broad_file1=None, broad_file2=None,
                      reorder='Diag', seabcmap='BuPu', sign=None, warn=None):

    """
    Saves and plots heatmap showing expression comparisons between all clusters of sp1 and all
    clusters of sp2.

    Args:
        result (numpy.array): matrix with weighted pearson correlation for gene expression
                              in all pairs of clusters
        cell_types1, cell_types2 (list): list of clusters in rows (species2) and columns (species2)
        outprefix (str): prefix for output files
        sp1, sp2 (str, optional): (short)name of species 1 and species 2 (to label axes on matrix)
        outformat (str, optional): format for saved figure (svg, png or pdf)
        broad_file1, broad_file2 (str, optional): path to pickle file with dict cluster to broad
                                                  annotations (file generated by pesci)
        reorder (str, optional): if not using broad annotation to order the heatmap, 3 alternatives
                                 are possible: 'Alpha' (alphabetical), 'Diag' will put best 
                                 1-1 match on the Diagonal, 'Clust' will use the order from
                                 hierarchical clustering using seaborn.clustermap
    """
    df = pd.DataFrame(data=result[0:,0:], columns=[sp2+'|'+i for i in cell_types2],
                      index=[sp1+'|'+i for i in cell_types1])
    df.to_csv(f'{outprefix}correlation_scores_matrix.tsv', sep='\t')

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

    # result[result<threshold_for_plot] = 0

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

    elif reorder == 'Diag':
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
        clustergrid = sns.clustermap(result, metric='correlation', method='average')

        row_idx_final = clustergrid.dendrogram_row.reordered_ind
        col_idx_final = clustergrid.dendrogram_col.reordered_ind

        result = result[row_idx_final,:]
        result = result[:,col_idx_final]
        cell_types2 = [list(cell_types2.keys())[i] for i in col_idx_final]
        cell_types1 = [list(cell_types1.keys())[i] for i in row_idx_final]
        plt.close('all')

    elif reorder != 'Alpha':
        logger.fatal('%s is not a valid argument for the reorder parameter, '
                       'use Diag, Clust or Alpha', reorder)
        raise ValueError("Argument Error")

    else:
        cell_types2 = list(cell_types2.keys())
        cell_types1 = list(cell_types1.keys())

    g = sns.clustermap(result, cmap=seabcmap, annot=False, vmin=0, vmax=1, xticklabels=cell_types2,
                yticklabels=cell_types1, cbar_kws={'label': 'weighted\ncorrelation', "shrink": 0.1},
                row_cluster=False, col_cluster=False, dendrogram_ratio=0.01,
                figsize=fsize, cbar_pos=(1.05, 0.8, 0.02, 0.15), row_colors=row_colors,
                col_colors=col_colors, colors_ratio=cratio)

    ax = g.ax_heatmap

    if sign is not None and not sign.empty:
        for pair in zip(sign[f"{sp1}_cell_cluster"], sign[f'{sp2}_cell_cluster']):
            idx_x = cell_types2.index(pair[1].split('|')[-1])
            idx_y = cell_types1.index(pair[0].split('|')[-1])
            ax.add_patch(Rectangle((idx_x, idx_y), 1, 1, fill=False, edgecolor='k', lw=1))

    if warn is not None and not warn.empty:
        for pair in zip(warn[f"{sp1}_cell_cluster"], warn[f'{sp2}_cell_cluster']):
            idx_x = cell_types2.index(pair[1].split('|')[-1])
            idx_y = cell_types1.index(pair[0].split('|')[-1])
            ax.text(idx_x + 0.5, idx_y + 0.6, '!', color='black',
                    path_effects=[pe.withStroke(linewidth=0.8, foreground="white")],
                    size=7, ha='center', va='center')

        ax.text(1.02, 0.7, '!: high score driven by\n< 10 genes', size=8,
                transform=plt.gcf().transFigure)

    if broad_file1 and broad_file2:
        handles = [Patch(facecolor=paldict[name]) for name in sorted(list(paldict.keys()))]
        plt.legend(handles, sorted(list(paldict.keys())), title='Broad annotation',
                   bbox_to_anchor=(1.26, 0.1), bbox_transform=plt.gcf().transFigure,
                   loc='lower right')

    ax.set_ylabel(sp1)
    ax.set_xlabel(sp2)
    plt.savefig(f'{outprefix}correlation_scores_matrix.{outformat}', bbox_inches='tight')
    plt.close('all')


def plot_expression_conservation(ec, n_ortho, outprefix):

    """
    Plots distribution of ec scores (histogram).

    Args:
        ec (numpy.array): vector with ec scores
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
                                 fc=1.5):

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
    coexp_number_enrich_score = []
    for i, clus1 in enumerate(mat1.clusters):
        for j, clus2 in enumerate(mat2.clusters):
            if result[i, j] > 0:
                genes_clus1_sp1 = np.where(mat1.matrix[:, i] > fc)[0]
                genes_clus2_sp2 = np.where(mat2.matrix[:, j] > fc)[0]
                ortho_idx = set(genes_clus1_sp1).intersection(set(genes_clus2_sp2))

                tot = len([k for k in ec if k>0])
                markers1 = len([k for k in genes_clus1_sp1 if ec[k]>0])
                markers2 = len([k for k in genes_clus2_sp2 if ec[k]>0])
                x = len([k for k in ortho_idx if ec[k]>0])

                # logger.info('%s-%s score = %s', sp1+'|'+clus1, sp2+'|'+clus2, result[i, j])
                # logger.info('Tot nb orthologs with ec>0: %s', tot)
                # logger.info('Nb Markers Sp1: %s', markers1)
                # logger.info('Nb Markers Sp2: %s', markers2)

                enr = 0
                if markers1 != 0 and markers2 != 0:
                    enr = x/(markers1*markers2/tot)

                coexp_number_enrich_score.append((sp1+'|'+clus1, sp2+'|'+clus2, x,
                                                  enr, result[i, j]))

                for k in ortho_idx:
                    g1, g2 = mat1.genes[k], mat2.genes[k]
                    fc1, fc2 = mat1.matrix[k, i], mat2.matrix[k, j]
                    ec_tmp = ec[k]
                    if ec_tmp > 0:
                        score_tmp = np.exp(np.log(fc1 + fc2)) * ec_tmp
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
    df.to_csv(f'{outprefix}gene_coexpression_table.tsv', sep='\t', index=False)

    df = pd.DataFrame.from_records(coexp_number_enrich_score,
                                   columns=[f'{sp1}_cell_cluster', f'{sp2}_cell_cluster',
                                            'n.co-expressed_genes', 'n.enrich',
                                            'weighted.correlation.score'])


    otsu_enrich = threshold_otsu(np.array([i[3] for i in coexp_number_enrich_score]))
    otsu_scores = threshold_otsu(np.array([i[4] for i in coexp_number_enrich_score]))

    logger.info('Automated threshold on weighted correlation scores: %s', round(otsu_scores, 3))

    logger.info('Automated threshold on co-expressed markers enrichment: %s', round(otsu_enrich, 3))


    sign = df.loc[(df['weighted.correlation.score'] > otsu_scores)
                   & (df['n.co-expressed_genes'] > 10) & (df['n.enrich'] > otsu_enrich)]
    sign = sign[[f'{sp1}_cell_cluster', f'{sp2}_cell_cluster']]

    warn = df.loc[(df['weighted.correlation.score'] > otsu_scores)
                   & (df['n.co-expressed_genes'] < 10)]

    warn = warn[[f'{sp1}_cell_cluster', f'{sp2}_cell_cluster']]

    df['automated.threshold.enrich'] = otsu_enrich
    df['automated.threshold.weighted.correlation.score'] = otsu_scores
    df.to_csv(f'{outprefix}coexpression_info.tsv', sep='\t', index=False)

    return warn, sign, otsu_scores, otsu_enrich





def compare(matrix_a, matrix_b, outprefix, sp1='sp1', sp2='sp2', random_id='',
            outformat='svg', min_fc=1.5, broad_file1=None,
            broad_file2=None, many_threshold=None, seabcmap='BuPu', use_thresh=False,
            plot_warn=True, reorder='Diag'):

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

    is_one2one = np.concatenate([np.ones(len(ortho)), np.zeros(len(para))])
    df = pd.DataFrame([genes1_ok, genes2_ok, ec, is_one2one]).T
    df.columns = ['genes sp1', 'genes sp2', 'expression conservation', 'is_one2one']
    df.to_csv(outprefix+random_id+sp1+'-'+sp2+'_'+'expression_conservation_scores.tsv',
              sep='\t', index=False)
    #compute weighted correlations between cell types of sp1 and cell types of sp2
    result = icc.column_wise_wcorr_einsum(mat1_ok, mat2_ok, ec)

    logger.info('Searching for co-expressed gene pairs')

    warn, sign, _, _ = make_coexpressed_genes_table(result, mat1, mat2, ec, len(ortho),
                                   outprefix+random_id+sp1+'-'+sp2+'_', sp1, sp2, fc=min_fc)

    if not plot_warn:
        warn = None

    if not use_thresh:
        sign = None

    logger.info('Saving outputs and plotting correlation matrix')

    plot_and_save_out(result, mat1.clusters, mat2.clusters,
                               outprefix+random_id+sp1+'-'+sp2+'_', sp1, sp2,
                               outformat=outformat,
                               broad_file1=broad_file1, broad_file2=broad_file2,
                               seabcmap=seabcmap, warn=warn, sign=sign, reorder=reorder)
