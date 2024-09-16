#!/usr/bin/env python3

"""
pesci compares gene expression of single-cell clusters between two species, using one-to-one and
many-to-many / many-to-one / one-to-many orthologes. pesci takes as input single-cell expression
count matrices and cluster annotations and uses the icc (iterative comparison of coexpression)
procedure to select orthologous gene pairs. Expression similarity between clusters is computed as a
weighted correlation.

Examples::

    $ pesci --matrix1 ../data/Cg_matrix_EM.tsv --matrix2 ../data/Pc_matrix_EM.tsv\
      --clusters1 ../data/Cragig_cell_id.tsv --clusters2 ../data/Procro_cell_id.tsv\
      --orthologous_gene_pairs ../data/orthologous_pairs.txt -c 10 -sp1 Cragig -sp2 Procro

    $ pesci --matrix1 ../data/placozoa/Tadh_final_matrix.tsv\
            --matrix2 ../data/placozoa/TrH2_final_matrix.tsv\
            --clusters1 ../data/placozoa/Tadh.sc_annot.tsv\
            --clusters2 ../data/placozoa/TrH2.sc_annot.tsv\
            -g ../data/placozoa/orthologous_pairs_ok.txt -c 10 -sp1 Tadh -sp2 TrH2
"""

import logging
import os
from pathlib import Path

import argparse

import coloredlogs

from . import normalize as nm, iterative_comparison_coexpression as icc, compare as cp


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 prog='pesci')

#Required arguments
required = parser.add_argument_group('Required arguments')
required.add_argument('-m1', '--matrix1', type=str, required=True,
                                        help="gene expression per cell for species 1 "
                                        "(raw count matrix), either: a dense matrix text file, "
                                        "a cellranger directory or an h5ad file")

required.add_argument('-m2', '--matrix2', type=str, required=True,
                                        help="gene expression per cell for species 2 "
                                        "(raw count matrix), either: a dense matrix text file, "
                                        "a cellranger directory or an h5ad file")

required.add_argument('-c1', '--clusters1', type=str, required=True,
                                          help="cell-to-clusters table for species 1 "
                                               "or name of cluster column in h5ad")

required.add_argument('-c2', '--clusters2', type=str, required=True,
                                          help="cell-to-clusters table for species 2 "
                                               "or name of cluster column in h5ad")

required.add_argument('-g', '--orthologous_gene_pairs', type=str, required=True,
                                                      help="gene orthologies file"
                                                      "(tab-delimited, 1 pair per line)")

#Optional arguments resources
resources = parser.add_argument_group('Resources')
resources.add_argument('-c', '--cores', help='Number of cores to use for parallel operations',
                                        type=int, required=False, default=1)


#Optional arguments running
ropt = parser.add_argument_group('Running options')

ropt.add_argument('--force', action='store_true',
                               help="recompute the per cluster normalized gene expression and ec "
                                    "scores even if output already exists")

ropt.add_argument('--seed', type=int, help="Fix random seed to ensure reproducible results",
                              default=123)

ropt.add_argument('--marker_specificity', type=float, help="Marker specificity (btwn 0.5 and 0.95)."
                                " Used in fold change expression calculation as "
                                "follows: 0.5 computes fold change in cluster 'a' vs median "
                                "expression across "
                                "clusters, 0.75 computes fold change in cluster 'a' vs 75 "
                                "expression percentile across clusters (3rd quartile)",
                                default=0.5)

ropt.add_argument('--ono2one_only', action='store_true', help="only use 1-to-1 orthologs")



ropt.add_argument('--ec_threshold_many', type=float, help="if provided, instead of best pair "
                                                "retain all distinct pairs with score > threshold",
                                            default=None)

ropt.add_argument('--random_id', type=str, help="triggers a randomization of orthologies and s"
                                                  "tore results in outdir/random_id", default='')

#Optional arguments input
iopt = parser.add_argument_group('Input options')

iopt.add_argument('--min_umi', type=int, help="Minimum total umi for a gene to be retained for "
                                                "comparison", default=10)

iopt.add_argument('--colclust', type=str, help="Name of the column corresponding to clusters in "
                                                 "the cell to cluster annotation file (does not "
                                                 "apply to h5ad input as it is specified by "
                                                 "--clusters1/--cluster2 directly). "
                                                 "If not provided, will use column `cluster_name`, "
                                                 "if it does not exist, will attempt to use the 2nd"
                                                 " column. Use --colclust1 and --colclust2 to "
                                                 "specify different column names for species1 and "
                                                 "species2.",
                    default='cluster_name')

iopt.add_argument('--colclust1', type=str, help="Name of the column corresponding to clusters in "
                                                 "the cell to cluster annotation file for species1 "
                                                 "(does not apply to h5ad input as it is specified "
                                                 "by --clusters1 directly). "
                                                 "If not provided, will use column `cluster_name`, "
                                                 "if it does not exist, will attempt to use the 2nd"
                                                 " column.",
                    default=None)

iopt.add_argument('--colclust2', type=str, help="Name of the column corresponding to clusters in "
                                                 "the cell to cluster annotation file for species2 "
                                                 "(does not apply to h5ad input as it is specified "
                                                 "by --cluster2 directly). "
                                                 "If not provided, will use column `cluster_name`, "
                                                 "if it does not exist, will attempt to use the 2nd"
                                                 " column.",
                    default=None)

iopt.add_argument('--colbroad', type=str, help="Name of the column corresponding to broad annot "
                                                 "in the cell to cluster annotation file or h5ad. "
                                                 "If not provided, will not use broad annotations "
                                                 "(this is for the heatmap plot only). Use "
                                                 "--colbroad1 and --colbroad2 to specify different "
                                                 "column names for species1 and species2.",
                    default=None)

iopt.add_argument('--colbroad1', type=str, help="Name of the column corresponding to broad annot "
                                                  "in the cell to cluster annotation file or h5ad "
                                                  "of species1. "
                                                  "If not provided, will not use broad annotations "
                                                  "(this is for the heatmap plot only).",
                    default=None)

iopt.add_argument('--colbroad2', type=str, help="Name of the column corresponding to broad annot "
                                                  "in the cell to cluster annotation file or h5ad "
                                                  "of species 2. "
                                                  "If not provided, will not use broad annotations "
                                                  "(this is for the heatmap plot only).",
                    default=None)

iopt.add_argument('--filter_out', type=str, required=False, default="",
                                    help='discard clusters starting with specified string,'
                                         'several can be specified as a comma-separated list '
                                         '(ex: --filter_out BC_,RGC_). Applies to both species1 '
                                         'and species2.')

iopt.add_argument('--filter_out1', type=str, required=False, default="",
                                    help='discard clusters starting with specified string,'
                                         'several can be specified as a comma-separated list '
                                         '(ex: --filter_out1 BC_,RGC_). Applies to species1 '
                                         'only.')

iopt.add_argument('--filter_out2', type=str, required=False, default="",
                                    help='discard clusters starting with specified string,'
                                         'several can be specified as a comma-separated list '
                                         '(ex: --filter_out2 BC_,RGC_). Applies to species2 '
                                         'only.')

iopt.add_argument('--keep_only', type=str, required=False, default="",
                                    help='keep only clusters starting with specified string, '
                                         'several can be specified as a comma-separated list '
                                         '(ex: --keep_only BC_,RGC_). Applies to both species1 '
                                         'and species2.')

iopt.add_argument('--keep_only1', type=str, required=False, default="",
                                    help='keep only clusters starting with specified string, '
                                         'several can be specified as a comma-separated list '
                                         '(ex: --keep_only1 BC_,RGC_). Applies to species1 '
                                         'only.')

iopt.add_argument('--keep_only2', type=str, required=False, default="",
                                    help='keep only clusters starting with specified string, '
                                         'several can be specified as a comma-separated list '
                                         '(ex: --keep_only2 BC_,RGC_)Applies to species2 '
                                         'only.')

#Optional arguments output
oopt = parser.add_argument_group('Output options')

oopt.add_argument('-o', '--outdir', type=str, required=False, default="output_pesci/",
                                    help='Name of output directory (will be created if does not '
                                         'exist)')

oopt.add_argument('-sp1', '--label_species1', type=str, required=False, default="sp1",
                                              help='Name of species1 (label for plots and outputs)')

oopt.add_argument('-sp2', '--label_species2', type=str, required=False, default="sp2",
                                              help='Name of species2 (label for plots and outputs)')

oopt.add_argument('-f', '--figure_format', type=str, required=False, default="svg",
                                           help='format for output figures',
                                           choices=['svg', 'png', 'pdf'])


oopt.add_argument('--plot_thresh', type=float, help="threshold to plot matches in resulting "
                                                      "comparison matrix & search for co-expressed"
                                                      " marker gene pairs", default=0)

oopt.add_argument('--min_fc', type=float, help="minimum fold-change to be considered marker of "
                                               "a cluster. Only used for the co-expressed marker "
                                               "gene table (to set in conjunction with "
                                               "--marker_specificity)", default=1.5)

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
norm_mat1 = args['outdir'] + 'files/' + sp1 + '_' + Path(MAT1).stem + '_expr_clusters_norm.tsv'
force_ec = False
if os.path.exists(norm_mat1) and os.path.getsize(norm_mat1) > 0 and not args['force']:
    logger.warning('Normalized gene-cluster expression matrix %s already exists '
                   'and will be used. Use --force to recompute.', norm_mat1)
else:
    colclust1 = args['colclust']
    if args['colclust1']:
        colclust1 = args['colclust1']
    broad1 = args['colbroad']
    if args['colbroad1']:
        broad1 =  args['colbroad1']

    keep_only = None
    if args['keep_only']:
        keep_only = tuple(args['keep_only'].split(','))
    if args['keep_only1']:
        keep_only = tuple(args['keep_only1'].split(','))

    filter_out_start = None
    if args['filter_out']:
        filter_out_start = tuple(args['filter_out'].split(','))
    if args['filter_out1']:
        filter_out_start = tuple(args['filter_out1'].split(','))

    nm.normalize(args['matrix1'], args['clusters1'], norm_mat1, cores=args['cores'],
                 filter_out_start=filter_out_start, keep_only=keep_only,
                 colclust=colclust1, broad=broad1, min_umi=args['min_umi'],
                 marker_specificity=args['marker_specificity'])
    force_ec = True

MAT2 = args['matrix2']
logger.info('Gene-cell expression matrix 2: %s', MAT2)
norm_mat2 = args['outdir'] + 'files/' + sp2 + '_' + Path(MAT2).stem + '_expr_clusters_norm.tsv'

if os.path.exists(norm_mat2) and os.path.getsize(norm_mat2) > 0 and not args['force']:
    logger.warning('Normalized gene-cluster expression matrix %s already exists'
                   ' and will be used. Use --force to recompute.', norm_mat2)
else:
    colclust2 = args['colclust']
    if args['colclust2']:
        colclust2 = args['colclust2']
    broad2 = args['colbroad']
    if args['colbroad2']:
        broad2 =  args['colbroad2']

    keep_only = None
    if args['keep_only']:
        keep_only = tuple(args['keep_only'].split(','))
    if args['keep_only2']:
        keep_only = tuple(args['keep_only2'].split(','))

    filter_out_start = None
    if args['filter_out']:
        filter_out_start = tuple(args['filter_out'].split(','))
    if args['filter_out2']:
        filter_out_start = tuple(args['filter_out2'].split(','))

    nm.normalize(args['matrix2'], args['clusters2'], norm_mat2, cores=args['cores'],
                 filter_out_start=filter_out_start, keep_only=keep_only,
                 colclust=colclust2, broad=broad2, min_umi=args['min_umi'],
                 marker_specificity=args['marker_specificity'])
    force_ec = True

if not args['random_id']:
    outprefix = args['outdir']+'files/'+sp1+'-'+sp2
    ec_scores = outprefix+'_1-to-1-orthologs_correlation_scores.tsv'

else:
    args['random_id'] = args['random_id'].rstrip('/') + '/'
    os.makedirs(args['outdir']+args['random_id'], exist_ok=True)
    outprefix = args['outdir']+args['random_id']+sp1+'-'+sp2
    ec_scores = outprefix+'_1-to-1-orthologs_correlation_scores.tsv'


if os.path.exists(ec_scores) and os.path.getsize(ec_scores) > 0 and not args['force']\
   and not force_ec:
    logger.warning('Orthologs expression conservation scores already computed %s'
                   ' and will be used. Homologs selection will not be re-run either.'
                   ' Use --force to recompute.', ec_scores)
else:
    icc.icc(norm_mat1, norm_mat2, args['orthologous_gene_pairs'], outprefix, max_combin=300,
            ncores=args['cores'], ono2one_only=args['ono2one_only'], random_id=args['random_id'],
            seed=args['seed'])

logger.info('Computing gene expression correlation between clusters of %s and %s', sp1, sp2)
broadfile1, broadfile2 = None, None
if args['colbroad'] or (args['colbroad1'] and args['colbroad2']):
    broadfile1 = norm_mat1.split('_matrix_')[0] + '_clusters_to_broad.pkl'
    broadfile2 = norm_mat2.split('_matrix_')[0] + '_clusters_to_broad.pkl'
cp.compare(norm_mat1, norm_mat2, args['outdir'], sp1, sp2, random_id=args['random_id'],
           threshold_for_plot=args['plot_thresh'], outformat=args['figure_format'],
           min_fc=args['min_fc'], broad_file1=broadfile1, broad_file2=broadfile2,
           many_threshold=args['ec_threshold_many'])
logger.info('Done! Results in %s', args['outdir'])
