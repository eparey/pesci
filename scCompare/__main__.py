#!/usr/bin/env python3

"""
python -m scCompare --matrix1 ../data/Cg_matrix_EM.csv --matrix2 ../data/Pc_matrix_EM.csv --clusters1 ../data/Cragig_cell_id.tsv --clusters2 ../data/Procro_cell_id.tsv --gene_families ../data/orthologous_pairs.txt -c 20 -sp1 Cragig -sp2 Procro

python -m scCompare --matrix1 ../data/placozoa/Tadh_final_matrix.tsv --matrix2 ../data/placozoa/TrH2_final_matrix.tsv --clusters1 ../data/placozoa/Tadh.sc_annot.tsv --clusters2 ../data/placozoa/TrH2.sc_annot.tsv --gene_families ../data/placozoa/orthologous_pairs_ok.txt -c 20 -sp1 Tadh -sp2 TrH2
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

parser.add_argument('-c1', '--clusters1', type=str, required=False, default= "leiden",
                                          help="cell-to-clusters table for species 1")

parser.add_argument('-c2', '--clusters2', type=str, required=False, default= "leiden",
                                          help="cell-to-clusters table for species 2")

parser.add_argument('-f', '--gene_families', type=str, required=True, help="gene families")

parser.add_argument('-c', '--cores', type=int, required=False, default=1)

parser.add_argument('-o', '--output_dir', type=str, required=False, default="out_scCompare/")

parser.add_argument('-sp1', '--label_species1', type=str, required=False, default="sp1")

parser.add_argument('-sp2', '--label_species2', type=str, required=False, default="sp2")

parser.add_argument('--filter_out', type=str, required=False, default="")

parser.add_argument('--force', action='store_true',
                               help="recompute the per cluster normalized gene expression")

parser.add_argument('--no_paralogs', action='store_true', help="only use 1-to-1 orthologs")

args = vars(parser.parse_args())

logger = logging.getLogger(__name__)
coloredlogs.install(level='INFO', logger=logger, fmt='%(asctime)s [%(levelname)s]: %(message)s',
                    field_styles={'levelname': {'color': ''}})

# Create output dir
args['output_dir'] = args['output_dir'].strip('/') + '/'
os.makedirs(args['output_dir']+ 'random', exist_ok=True)

sp1 = args["label_species1"]
sp2 = args["label_species2"]

MAT1 = args['matrix1']
logger.info(f'Gene-cell expression matrix 1: {MAT1}')
norm_mat1 = args['output_dir'] + Path(MAT1).stem + '_expr_clusters_norm.tsv'

if os.path.exists(norm_mat1) and os.path.getsize(norm_mat1) > 0 and not args['force']:
    logger.warning(f'Normalized gene-cluster expression matrix {norm_mat1} already exists\
                     and will be used. Use --force to recompute.')
else:
    nm.normalize_main(args['matrix1'], args['clusters1'], norm_mat1, cores=args['cores'],
                      filter_out_start=args['filter_out'])

MAT2 = args['matrix2']
logger.info(f'Gene-cell expression matrix 2: {MAT2}')
norm_mat2 = args['output_dir'] + Path(MAT2).stem + '_expr_clusters_norm.tsv'

if os.path.exists(norm_mat2) and os.path.getsize(norm_mat2) > 0 and not args['force']:
    logger.warning(f'Normalized gene-cluster expression matrix {norm_mat2} already exists\
                     and will be used. Use --force to recompute.')
else:
    nm.normalize_main(args['matrix2'], args['clusters2'], norm_mat2,
                      filter_out_start=args['filter_out'])

icc.icc_main(norm_mat1, norm_mat2, args['gene_families'],
                        args['output_dir']+sp1+'-'+sp2+'_', max_combin=200, ncores=args['cores'],
                        noparalogs=args['no_paralogs'])

logger.info(f'Computing gene expression correlation between clusters of {sp1} and {sp2}')
cp.compare_main(norm_mat1, norm_mat2, args['output_dir']+sp1+'-'+sp2+'_', args['cores'],
                args['label_species1'], args['label_species2'])
OUTDIR = args['output_dir']
logger.info(f'Done! Results in {OUTDIR}')
