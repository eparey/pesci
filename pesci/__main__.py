#!/usr/bin/env python3

"""
pesci command-line entry point.

Takes as input single-cell expression count matrices and cluster annotations for 2 species and uses
the icc (iterative correlation coexpression) procedure to select orthologous gene pairs and compute
expression similarity between clusters as a weighted correlation.

Examples::

    $ pesci --matrix1 ../data/Cg_matrix_EM.tsv --matrix2 ../data/Pc_matrix_EM.tsv \
      --clusters1 ../data/Cragig_cell_id.tsv --clusters2 ../data/Procro_cell_id.tsv \
      --orthologous_gene_pairs ../data/orthologous_pairs.txt -c 10 -sp1 Cragig -sp2 Procro

    $ pesci --matrix1 ../data/placozoa/Tadh_final_matrix.tsv \
            --matrix2 ../data/placozoa/TrH2_final_matrix.tsv \
            --clusters1 ../data/placozoa/Tadh.sc_annot.tsv \
            --clusters2 ../data/placozoa/TrH2.sc_annot.tsv \
            -g../data/placozoa/orthologous_pairs_ok.txt -c 10 -sp1 Tadh -sp2 TrH2
"""

import logging
import os
from pathlib import Path

import argparse

import coloredlogs

from . import normalize as nm, iterative_comparison_coexpression as icc, compare as cp


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-m1', '--matrix1', type=str, required=True,
                                        help="gene expression per cell for species 1")

parser.add_argument('-m2', '--matrix2', type=str, required=True,
                                        help="gene expression per cell for species 1")

parser.add_argument('-c1', '--clusters1', type=str, required=True,
                                          help="cell-to-clusters table for species 1 "
                                               "or name of cluster column in h5ad")

parser.add_argument('-c2', '--clusters2', type=str, required=True,
                                          help="cell-to-clusters table for species 2 "
                                               "or name of cluster column in h5ad")

parser.add_argument('-g', '--orthologous_gene_pairs', type=str, required=True,
                                                      help="gene orthologies "
                                                      "(tab-delimited, 1 pair per line)")

parser.add_argument('-c', '--cores', type=int, required=False, default=1)

parser.add_argument('-o', '--outdir', type=str, required=False, default="output_pesci/")

parser.add_argument('-sp1', '--label_species1', type=str, required=False, default="sp1")

parser.add_argument('-sp2', '--label_species2', type=str, required=False, default="sp2")

parser.add_argument('-f', '--figure_format', type=str, required=False, default="svg",
                                             choices=['svg', 'png', 'pdf'])

parser.add_argument('--filter_out', type=str, required=False, default="",
                                    help='discard clusters starting with specified string')

parser.add_argument('--force', action='store_true',
                               help="recompute the per cluster normalized gene expression and ec "
                                    "scores")

parser.add_argument('--ono2one_only', action='store_true', help="only use 1-to-1 orthologs")

parser.add_argument('--random_id', type=str, help="triggers a randomization of orthologies and s"
                                                  "tore results in outdir/random_id", default='')

parser.add_argument('--plot_thresh', type=float, help="threshold to plot matches in resulting "
                                                      "comparison matrix & search for co-expressed"
                                                      " marker gene pairs", default=0)

parser.add_argument('--min_fc', type=float, help="minimum fold-change ot be considered marker of "
                                                 "a cluster. Only used for co-expressed marker "
                                                 "gene table.", default=1.5)

parser.add_argument('--seed', type=int, help="Fix random seed to ensure reproducible results",
                              default=123)


args = vars(parser.parse_args())

logger = logging.getLogger(__name__)
coloredlogs.install(level='INFO', logger=logger, fmt='%(asctime)s [%(levelname)s]: %(message)s',
                    field_styles={'levelname': {'color': ''}})

# Create output dir
args['outdir'] = args['outdir'].rstrip('/') + '/'
os.makedirs(args['outdir']+ 'files', exist_ok=True)

sp1 = args["label_species1"]
sp2 = args["label_species2"]

MAT1 = args['matrix1']
logger.info('Gene-cell expression matrix 1: %s', MAT1)
norm_mat1 = args['outdir'] + 'files/' + sp2 + '_' + Path(MAT1).stem + '_expr_clusters_norm.tsv'

if os.path.exists(norm_mat1) and os.path.getsize(norm_mat1) > 0 and not args['force']:
    logger.warning('Normalized gene-cluster expression matrix %s already exists '
                   'and will be used. Use --force to recompute.', norm_mat1)
else:
    nm.normalize(args['matrix1'], args['clusters1'], norm_mat1, cores=args['cores'],
                 filter_out_start=args['filter_out'])

MAT2 = args['matrix2']
logger.info('Gene-cell expression matrix 2: %s', MAT2)
norm_mat2 = args['outdir'] + 'files/' + sp2 + '_' + Path(MAT2).stem + '_expr_clusters_norm.tsv'

if os.path.exists(norm_mat2) and os.path.getsize(norm_mat2) > 0 and not args['force']:
    logger.warning('Normalized gene-cluster expression matrix %s already exists'
                   ' and will be used. Use --force to recompute.', norm_mat2)
else:
    nm.normalize(args['matrix2'], args['clusters2'], norm_mat2,
                      filter_out_start=args['filter_out'])

if not args['random_id']:
    outprefix = args['outdir']+'files/'+sp1+'-'+sp2
    ec_scores = outprefix+'_1-to-1-orthologs_correlation_scores.tsv'

else:
    args['random_id'] = args['random_id'].rstrip('/') + '/'
    os.makedirs(args['outdir']+args['random_id'], exist_ok=True)
    outprefix = args['outdir']+args['random_id']+sp1+'-'+sp2
    ec_scores = outprefix+'_1-to-1-orthologs_correlation_scores.tsv'


if os.path.exists(ec_scores) and os.path.getsize(ec_scores) > 0 and not args['force']:
    logger.warning('Orthologs expression conservation scores already computed %s'
                   ' and will be used. Homologs selection will not be re-run either.'
                   ' Use --force to recompute.', ec_scores)
else:
    icc.icc(norm_mat1, norm_mat2, args['gene_families'], outprefix, max_combin=300,
            ncores=args['cores'], ono2one_only=args['ono2one_only'], random_id=args['random_id'],
            seed=args['seed'])

logger.info('Computing gene expression correlation between clusters of %s and %s', sp1, sp2)
cp.compare(norm_mat1, norm_mat2, args['outdir'], sp1, sp2, random_id=args['random_id'],
           threshold_for_plot=args['plot_thresh'], outformat=args['figure_format'],
           min_fc=args['min_fc'])
logger.info('Done! Results in %s', args['outdir'])
