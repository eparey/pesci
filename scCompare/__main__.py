import sys

from pathlib import Path
import os
import argparse

import coloredlogs, logging

from . import normalize as nm, select_paralogs as sp, compare as cp


if __name__ == '__main__':
    

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-m1', '--matrix1', type=str, required=True, help="gene expression per cell for species 1")

    PARSER.add_argument('-m2', '--matrix2', type=str, required=True, help="gene expression per cell for species 1")

    PARSER.add_argument('-c1', '--clusters1', type=str, required=True, help="cell-to-clusters table for species 1")

    PARSER.add_argument('-c2', '--clusters2', type=str, required=True, help="cell-to-clusters table for species 2")

    PARSER.add_argument('-f', '--gene_families', type=str, required=True, help="gene families")

    PARSER.add_argument('-c', '--cores', type=int, required=False, default=1)

    PARSER.add_argument('-o', '--output_dir', type=str, required=False, default="out_scCompare/")

    PARSER.add_argument('-sp1', '--label_species1', type=str, required=False, default="sp1")

    PARSER.add_argument('-sp2', '--label_species2', type=str, required=False, default="sp2")

    PARSER.add_argument('--filter_out', type=str, required=False, default="")

    PARSER.add_argument('--random', action='store_true', help="this triggers a randomization of the gene homologies to evaluate the expected random background correlations")

    # PARSER.add_argument('--reoptimize', action='store_true', help="this triggers a re-optimization of ec scores on the whole set after paralogs selection")

    PARSER.add_argument('--no_paralogs', action='store_true', help="only use 1-to-1 orthologs")

    ARGS = vars(PARSER.parse_args())

    logger = logging.getLogger(__name__)
    coloredlogs.install(level='INFO', logger=logger, fmt='%(asctime)s - [%(levelname)s]: %(message)s', field_styles={'levelname': {'color': ''}})

    # Create output dir
    ARGS['output_dir'] = ARGS['output_dir'].strip('/') + '/' 
    os.makedirs(ARGS['output_dir']+ 'random', exist_ok=True)

    sp1 = ARGS["label_species1"]
    sp2 = ARGS["label_species2"]

    MAT1 = ARGS['matrix1']
    logger.info(f'Gene-cell expression matrix 1: {MAT1}')
    norm_mat1 = ARGS['output_dir'] + Path(MAT1).stem + '_expr_clusters_norm.tsv'

    if os.path.exists(norm_mat1) and os.path.getsize(norm_mat1) > 0:
        logger.warning(f'Normalized gene-cluster expression matrix {norm_mat1} already exists and will be used. Use --force to recompute.')
    else:
        nm.normalize_main(ARGS['matrix1'], ARGS['clusters1'], norm_mat1, cores=ARGS['cores'], filter_out_start=ARGS['filter_out'])

    MAT2 = ARGS['matrix2']
    logger.info(f'Gene-cell expression matrix 2: {MAT2}')
    norm_mat2 = ARGS['output_dir'] + Path(MAT2).stem + '_expr_clusters_norm.tsv'

    if os.path.exists(norm_mat2) and os.path.getsize(norm_mat2) > 0:
        logger.warning(f'Normalized gene-cluster expression matrix {norm_mat2} already exists and will be used. Use --force to recompute.')
    else:
        nm.normalize_main(ARGS['matrix2'], ARGS['clusters2'], norm_mat2, filter_out_start=ARGS['filter_out'])

    sp.select_paralogs_main(norm_mat1, norm_mat2, ARGS['gene_families'], ARGS['output_dir']+sp1+'-'+sp2+'_', max_combin=200,
                            ncores=ARGS['cores'], noparalogs=ARGS['no_paralogs'])

    logger.info(f'Computing gene expression correlation between clusters of {sp1} and clusters of {sp2}')
    cp.compare_main(norm_mat1, norm_mat2, ARGS['output_dir']+sp1+'-'+sp2+'_', ARGS['cores'], ARGS['label_species1'], ARGS['label_species2']) #, random=random_homologs_dir reoptimize_ec=ARGS['reoptimize'] #remove the reoptimize arg
    OUTDIR = ARGS['output_dir']
    logger.info(f'Done! Results in {OUTDIR}')

"""
python -m scCompare --matrix1 ../data/Cg_matrix_EM.csv --matrix2 ../data/Pc_matrix_EM.csv --clusters1 ../data/Cragig_cell_id.tsv --clusters2 ../data/Procro_cell_id.tsv --gene_families ../data/orthologous_pairs.txt -c 20 -sp1 Cragig -sp2 Procro
"""

"""
python -m scCompare --matrix1 ../data/placozoa/Tadh_final_matrix.tsv --matrix2 ../data/placozoa/TrH2_final_matrix.tsv --clusters1 ../data/placozoa/Tadh.sc_annot.tsv --clusters2 ../data/placozoa/TrH2.sc_annot.tsv --gene_families ../data/placozoa/orthologous_pairs_ok.txt -c 20 -sp1 Tadh -sp2 TrH2
"""