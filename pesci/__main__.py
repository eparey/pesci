#!/usr/bin/env python3

"""
Pesci compares gene expression of single-cell clusters between two species using the ICC algorithm,
using one-to-one orthologues complemented with many-to-many / many-to-one / one-to-many orthologous
genes.

Examples::

    $ pesci -m1 data/Cragig_matrix_EM.tsv.gz -m2 data/Procro_matrix_EM.tsv.gz\
      -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv\
      -g data/orthologous_pairs_Procro-Cragig.txt -c 10 -sp1 Cragig -sp2 Procro

    $ pesci --matrix1 ../data/placozoa/Tadh_final_matrix.tsv.gz\
            --matrix2 ../data/placozoa/TrH2_final_matrix.tsv.gz\
            --clusters1 ../data/placozoa/Tadh.sc_annot.tsv\
            --clusters2 ../data/placozoa/TrH2.sc_annot.tsv\
            -g ../data/placozoa/orthologous_pairs_ok.txt -c 10 -sp1 Tadh -sp2 TrH2

"""

import sys
import logging
import os
from pathlib import Path
from functools import partialmethod

import argparse

import traceback
import coloredlogs


import tqdm

from . import version
from . import normalize as nm, iterative_comparison_coexpression as icc, compare as cp

#do not show progress bar if stderr is redirected to file
if not sys.stderr.isatty():
    tqdm.tqdm.__init__ = partialmethod(tqdm.tqdm.__init__, disable=True)

def parse_commandline():
    """
    pesci command-line argument parser
    """
    parser = argparse.ArgumentParser(description="Pesci compares gene expression of single-cell "
                                     "clusters between two species using the ICC algorithm, using "
                                     "one-to-one orthologues complemented with many-to-many / "
                                     "many-to-one / one-to-many orthologous genes.\n\n"
                                     "Example usage:\npesci -m1 data/Cg_matrix_EM.tsv "
                                     "-m2 data/Pc_matrix_EM.tsv -c1 data/Cragig_cell_id.tsv "
                                     "-c2 data/Procro_cell_id.tsv "
                                     "-g data/orthologous_pairs_Procro-Cragig.txt -c 10 "
                                     "-sp1 Cragig -sp2 Procro",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     prog='pesci')

    parser.add_argument('-v', '--version', action='version',
                                           version='%(prog)s ' + version.__version__)

    #Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-m1', '--matrix1', type=str, nargs='+', required=True,
                                              help="gene expression per cell for species 1 "
                                              "(raw count matrix), accepts: dense matrix (can "
                                              "be .gz), sparse matrix (CellRanger-like), or "
                                              "scanpy `.h5ad`. Multiple files allowed if "
                                              "several libraries.")

    required.add_argument('-m2', '--matrix2', type=str, nargs='+', required=True,
                                              help="gene expression per cell for species 2, "
                                                   "same format options as --matrix1")

    required.add_argument('-c1', '--clusters1', type=str, required=True,
                                              help="cell-to-clusters table for species 1 "
                                                   "or name of cluster column in h5ad")

    required.add_argument('-c2', '--clusters2', type=str, required=True,
                                              help="cell-to-clusters table for species 2 "
                                                   "or name of cluster column in h5ad")

    required.add_argument('-g', '--ortho_pairs', type=str, required=False,
                                                           help="gene orthology"
                                                           " file (one pair per line).")

    #Recommended optional arguments
    raopt = parser.add_argument_group('recommended arguments')
    raopt.add_argument('-o', '--outdir', type=str, required=False, default="output_pesci/",
                                        help='name of output directory (will be created if does '
                                             'not exist)')

    raopt.add_argument('-l1', '--label1', type=str, required=False, default="sp1",
                                                  help='name of species2/dataset1 (label for '
                                                       'plots & outputs)')


    raopt.add_argument('-l2', '--label2', type=str, required=False, default="sp2",
                                                  help='name of species2/dataset2 (label for '
                                                  'plots & outputs)')

    raopt.add_argument('--within_species', action='store_true',
                                           help="use this flag when comparing datasets "
                                                "within the same species")

    #Optional arguments resources
    resources = parser.add_argument_group('resources')
    resources.add_argument('-c', '--cores', help='number of cores to use for parallel operations',
                                            type=int, required=False, default=1)


    #Optional arguments running
    ropt = parser.add_argument_group('running')

    ropt.add_argument('--force', action='store_true',
                                   help="recompute all intermediary files (per cluster normalized "
                                   "gene expression and ec scores) even if outputs already exist")

    ropt.add_argument('--seed', type=int, help="set random seed", default=123)

    ropt.add_argument('--marker_specificity', type=float, help="marker specificity "
                                    "(btwn 0.5 and 0.95) - used in fold change expression "
                                    "calculation as follows: 0.5 computes fold change in cluster "
                                    "'a' vs median expression across clusters, 0.75 computes fold "
                                    "change in cluster 'a' vs 75 expression percentile across "
                                    "clusters", default=0.5)

    ropt.add_argument('--ono2one_only', action='store_true', help="use 1-to-1 orthologs only")

    ropt.add_argument('--ec_threshold_many', type=float, help="if set, instead of using only "
                                             "the best pair for many-to-many and "
                                             "many-to-one orthologs, "
                                             "retain all distinct pairs with score > threshold",
                                             default=None)

    ropt.add_argument('--do_not_downsample', action='store_true', help="use all 1-to-1 orthologs "
                           "(instead of 1000 randomly chosen) for selection of best pairs for "
                           " the many-to-many and many-to-one orthologs (slower)")

    ropt.add_argument('--random_id', type=str, help="triggers a randomization of orthologies "
                                                    "to run pesci on randomized orthologs; "
                                                    "stores results in outdir/random_id", 
                                                    default='')

    ropt.add_argument('--max_combin', type=int, help="Maximum combination for many-to-many "
                                                     "ortholog groups, highly mutligenic "
                                                     "family will be skipped",
                                                     default=300)

    ropt.add_argument('--no_pbar', action='store_true',
                                   help="do not show progress bar")

    #Optional arguments input
    iopt = parser.add_argument_group('inputs')

    iopt.add_argument('--min_umi', type=int, help="minimum total umi for a gene to be retained for "
                                                  "comparison", default=10)

    iopt.add_argument('--colclust', type=str, help="name of the column corresponding to clusters "
                                                   "in the cell to cluster annotation file (does "
                                                   "not apply to h5ad input as it is specified by "
                                                   "--clusters1/--cluster2 directly) - "
                                                   "if not set, will use column `cluster_name`, "
                                                   "if it does not exist, will attempt to use the "
                                                   "2nd column \n"
                                                   "use --colclust1 or --colclust2 to "
                                                   "specify different column names for species1 or"
                                                   " species2.",
                                    default='cluster_name')

    iopt.add_argument('--colclust1', type=str, help="name of the column corresponding to clusters "
                                                    "in the cell to cluster annotation file for "
                                                    "species1 (does not apply to h5ad input as it "
                                                    "is specified by --clusters1 directly) - "
                                                    "if not set, will use column `cluster_name`, "
                                                    "if it does not exist, will attempt to use the "
                                                    "2nd column.",
                                     default=None)

    iopt.add_argument('--colclust2', type=str, help="name of the column corresponding to clusters "
                                                    "in the cell to cluster annotation file for "
                                                    "species2 (does not apply to h5ad input as it "
                                                    "is specified by --clusters1 directly) - "
                                                    "if not set, will use column `cluster_name`, "
                                                    "if it does not exist, will attempt to use the "
                                                    "2nd column.",
                                     default=None)

    iopt.add_argument('--colbroad', type=str, help="name of the column corresponding to broad "
                                                    "annotations in the cell to cluster annotation "
                                                    "file or h5ad - if not set, will not use "
                                                    "broad annotations (this is for the heatmap "
                                                    "plot only) \n"
                                                    "use --colbroad1 or --colbroad2 "
                                                    "to specify different column names for species1"
                                                    " or species2.",
                                    default=None)

    iopt.add_argument('--colbroad1', type=str, help="name of the column corresponding to broad "
                                                    "annotations in the cell to cluster annotation "
                                                    "file or h5ad of species1 - if not set, "
                                                    "will not use broad annotations (this is for "
                                                    "the heatmap plot only).",
                                     default=None)

    iopt.add_argument('--colbroad2', type=str, help="name of the column corresponding to broad "
                                                    "annotations in the cell to cluster annotation "
                                                    "file or h5ad of species2 - if not set, "
                                                    "will not use broad annotations (this is for "
                                                    "the heatmap plot only).",
                                     default=None)

    iopt.add_argument('--filter_out', type=str, required=False, default="",
                                      help='discard clusters starting with specified string,'
                                           'several can be specified as a comma-separated list '
                                           '(ex: --filter_out BC_,RGC_). Applies to both species1 '
                                           '& species2.')

    iopt.add_argument('--filter_out1', type=str, required=False, default="",
                                        help='discard clusters starting with specified string,'
                                             'several can be specified as a comma-separated list '
                                             '(ex: --filter_out1 BC_,RGC_) - applies to species1 '
                                             'only.')

    iopt.add_argument('--filter_out2', type=str, required=False, default="",
                                        help='discard clusters starting with specified string,'
                                             'several can be specified as a comma-separated list '
                                             '(ex: --filter_out2 BC_,RGC_) - applies to species2 '
                                             'only.')

    iopt.add_argument('--keep_only', type=str, required=False, default="",
                                        help='keep only clusters starting with specified string, '
                                             'several can be specified as a comma-separated list '
                                             '(ex: --keep_only BC_,RGC_) - applies to both species1'
                                             ' & species2.')

    iopt.add_argument('--keep_only1', type=str, required=False, default="",
                                        help='keep only clusters starting with specified string, '
                                             'several can be specified as a comma-separated list '
                                             '(ex: --keep_only1 BC_,RGC_) - applies to species1 '
                                             'only.')

    iopt.add_argument('--keep_only2', type=str, required=False, default="",
                                        help='keep only clusters starting with specified string, '
                                             'several can be specified as a comma-separated list '
                                             '(ex: --keep_only2 BC_,RGC_) - applies to species2 '
                                             'only.')

    #Optional arguments output
    oopt = parser.add_argument_group('outputs')

    oopt.add_argument('-l', '--logfile', type=str, required=False, default="",
                                               help='logfile, if specified this will overwrite'
                                               ' the default (i.e. pesci.log in the specified'
                                               ' output folder)')

    oopt.add_argument('-f', '--figure_format', type=str, required=False, default="pdf",
                                               help='format for output figures',
                                               choices=['svg', 'png', 'pdf'])

    oopt.add_argument('-r', '--reorder', type=str, required=False, default="Diag",
                                                  help='method for ordering cell clusters on the '
                                                  'heatmap: Diag (default): try to maximise '
                                                  ' matches on the diagonal, keeping input order '
                                                  'as much as possible; '
                                                  'Clust: use hierarchical clustering of rows and '
                                                  'columns (average linkage, euclidean distance) '
                                                  'or Alpha: alphabetical. This option '
                                                  'is ignored if providing broad annotations.',
                                                  choices=['Diag','Clust', 'Alpha'])

    oopt.add_argument('--min_fc', type=float, help="minimum fold-change to be considered marker of "
                                                   "a cluster - only used for the co-expressed "
                                                   "marker genes table (to set in conjunction with "
                                                   "--marker_specificity).", default=1.5)


    oopt.add_argument('--seaborn_cmap', type=str, required=False, default="BuPu",
                                        help='name of the seaborn colormap for the heatmap.')

    oopt.add_argument('--show_auto_threshold', action='store_true', help="experimental option to "
                                                                        "highlight matches above "
                                                                        "thresholds on the "
                                                                        "heatmap plot.")
    oopt.add_argument('--do_not_plot_warn', action='store_false', help="do not plot warning on "
                                                                        "the heatmap for matches "
                                                                        "driven by < 10 genes.")

    args = vars(parser.parse_args())

    return args


def configure_logs(outdir, logfile=None):
    """
    pesci logs
    """
    logger = logging.getLogger()
    if not logfile:
        logfile = outdir + 'pesci.log'
    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    def log_uncaught_exc(etype, evalue, tb):
        logging.error(''.join(traceback.format_tb(tb)))
        if str(evalue):
            logging.error('%s: %s', etype.__name__, evalue)
        else:
            logging.error('%s', etype.__name__)

    sys.excepthook = log_uncaught_exc

    with open(outdir + 'pesci.log', 'a', encoding="utf-8") as outlogfile:
        outlogfile.write('# '+" ".join(sys.argv)+'\n')
        outlogfile.write('# '+'pesci ' + version.__version__+'\n')

    coloredlogs.install(level='INFO', logger=logger, fmt='%(asctime)s [%(levelname)s]: %(message)s',
                        field_styles={'levelname': {'color': ''}})

    return logger


def main():
    """
    pesci cli main
    """

    args = parse_commandline()

    if args['no_pbar']:
        tqdm.tqdm.__init__ = partialmethod(tqdm.tqdm.__init__, disable=True)

    # Create output dir
    args['outdir'] = args['outdir'].rstrip('/') + '/'
    os.makedirs(args['outdir']+ 'files', exist_ok=True)
    outdir = args['outdir']

    logger = configure_logs(outdir, args['logfile'])

    if not (args['ortho_pairs'] or args['within_species']):
        logger.fatal('Either provide an orthology file with --ortho_pairs or '
                     'use --within_species for comparison across organs of the same species')

    sp1 = args["label1"]
    sp2 = args["label2"]

    mat1 = args['matrix1']
    logger.info('Gene-cell expression matrix 1: %s', ' '.join(mat1))
    norm_mat1 = args['outdir'] + 'files/' + sp1 + '_' + Path(mat1[0]).stem +\
               '_expr_clusters_norm_all.tsv'
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

    mat2 = args['matrix2']
    logger.info('Gene-cell expression matrix 2: %s', ' '.join(mat2))
    norm_mat2 = args['outdir'] + 'files/' + sp2 + '_' + Path(mat2[0]).stem +\
               '_expr_clusters_norm_all.tsv'

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
        ec_scores_para = outprefix+'_orthologs_many_correlation_scores.txt'

    else:
        args['random_id'] = args['random_id'].rstrip('/') + '/'
        os.makedirs(args['outdir']+args['random_id'], exist_ok=True)
        outprefix = args['outdir']+args['random_id']+sp1+'-'+sp2
        ec_scores = outprefix+'_1-to-1-orthologs_correlation_scores.tsv'
        ec_scores_para = outprefix+'_orthologs_many_correlation_scores.txt'
        args['seed'] = None


    if os.path.exists(ec_scores) and os.path.getsize(ec_scores) > 0 and not args['force']\
       and not force_ec and\
       ((os.path.exists(ec_scores_para) and os.path.getsize(ec_scores_para) > 0)\
        or args['ono2one_only']):
        logger.warning('Orthologs expression conservation scores already computed %s'
                       ' and will be used. Homologs selection will not be re-run either.'
                       ' Use --force to recompute.', ec_scores)
    else:
        icc.icc(norm_mat1, norm_mat2, args['ortho_pairs'], outprefix,
                max_combin=args['max_combin'],
                ncores=args['cores'], ono2one_only=args['ono2one_only'],
                random_id=args['random_id'], seed=args['seed'],
                do_not_downsample=args['do_not_downsample'],
                within_species=args['within_species'])

    logger.info('Computing gene expression correlation between clusters of %s and %s', sp1, sp2)
    broadfile1, broadfile2 = None, None
    if args['colbroad'] or (args['colbroad1'] and args['colbroad2']):
        broadfile1 = norm_mat1.split('_matrix_')[0] + '_clusters_to_broad.pkl'
        broadfile2 = norm_mat2.split('_matrix_')[0] + '_clusters_to_broad.pkl'

    cp.compare(norm_mat1, norm_mat2, args['outdir'], sp1, sp2, random_id=args['random_id'],
               outformat=args['figure_format'], min_fc=args['min_fc'], broad_file1=broadfile1,
               broad_file2=broadfile2, many_threshold=args['ec_threshold_many'],
               seabcmap=args['seaborn_cmap'], use_thresh=args['show_auto_threshold'],
               plot_warn=args['do_not_plot_warn'], reorder=args['reorder'])

    logger.info('Done! Results in %s', args['outdir'])


if __name__ == '__main__':
    main()
