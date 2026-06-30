"""
Pesci utility script to summarize pairwise cross-species comparisons stored in a
Pesci output folder. Uses hierarchical clustering (avg linkage) on1 - correlation distance matrix to visualise distances between cell clusters.
"""

import logging

import glob

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import scipy.cluster.hierarchy as sch

import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install()


def list_files(pesci_outdir):

    """
    List all Pesci correlation matrix files present in supplied directory.
    Args:
        pesci_outdir (str): pesci directory

    Returns:
        list of pd.DataFrame: list of correlation matrices
    """

    infiles = glob.glob(f'{pesci_outdir}/*_correlation_scores_matrix.tsv')
    species = {tuple(i.split('/')[-1].split('_correlation_scores_matrix.tsv')[0].split('-')) 
               for i in infiles}
    species = {element for item in species for element in item}

    logger.info('Found %s pairwise comparison files in %s, '
                'for a total of %s different datasets (%s).',
                len(infiles), pesci_outdir, len(species), ', '.join(species))

    expected = len(species) * (len(species) + 1) / 2
    if len(infiles) != expected:
        logger.fatal('Expecting exactly %s files for all pairwise comparisons '
                     '(including self-comparisons), but found %s (%s).', 
             int(expected), len(infiles), ', '.join(infiles))
        raise ValueError("Input Error")


    df_list = []
    for i in infiles:
        df = pd.read_csv(i, sep='\t', index_col=0)

        df = df.clip(lower=0, upper=1.0)

        df_list.append(df)

    return df_list




def combine_sim_mat(df_list):

    """
    Combines a list of pesci correlation matrix (pd.Dataframe) into a single matrix.
    Args:
        df_list (list of pd.DataFrame): data to combine

    Returns:
        pd.DataFrame: all vs all correlation matrix (i.e. all pairwise comparisons)
    """

    long_forms = []

    for df in df_list:
        # Convert matrix to a long-format table: [Row, Column, Distance]
        melted = df.melt(ignore_index=False).reset_index()
        melted.columns = ["cell_cluster", "cell_cluster_col", "dist"]

        long_forms.append(melted)

        if set(df.index) != set(df.columns):
            inverse = melted.rename(columns={"cell_cluster": "cell_cluster_col", 
                                             "cell_cluster_col": "cell_cluster"})
            long_forms.append(inverse)

    master_long = pd.concat(long_forms, ignore_index=True)

    final_matrix = master_long.pivot(
        index="cell_cluster", columns="cell_cluster_col", values="dist"
    )

    final_matrix = final_matrix.reindex(
        index=final_matrix.index, columns=final_matrix.index
    )

    return final_matrix


def hierarchical_clustering(pesci_outdir, jobname, outformat='pdf'):

    """
    Transform correlation matrix to a distance matrix (1 - corr) and performs hierarchical 
    clustering. Output figure and correlation matrix are stored to files in the supplied directory.

    Args:
        pesci_outdir (str): pesci directory
        jobname (str): prefix for the output files
        outformat (str, optional): format for saved figure (svg, png or pdf)

    """

    df_list = list_files(pesci_outdir)

    df = combine_sim_mat(df_list)

    df.to_csv(f'{pesci_outdir}/{jobname}_all_vs_all_pesci_scores.tsv', sep='\t')

    df = 1.0 - df

    dist_matrix = np.array(df.values, copy=True)
    condensed_distances = sch.distance.squareform(dist_matrix, checks=False)

    clus = sch.linkage(condensed_distances, method="average")

    n = len(clus) + 1
    height = max(6, min(8, n * 0.05 + 3.5))
    _, ax = plt.subplots(figsize=(8, height))

    with plt.rc_context({"lines.linewidth": 0.8}):
        sch.dendrogram(
        clus,
        labels=df.index.tolist(),
        orientation="left", 
        leaf_font_size=6,
        color_threshold=0,
        above_threshold_color="black",
        ax=ax,
        )


    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.get_xaxis().set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{pesci_outdir}/{jobname}_all_vs_all_hierarchical_clustering.{outformat}')

    plt.close('all')
