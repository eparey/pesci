from pathlib import Path
import os
import argparse

from . import normalize as nm, select_paralogs as sp, compare as cp


PARSER = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

PARSER.add_argument('-m1', '--matrix1', type=str, required=True, help="gene expression per cell for species 1")

PARSER.add_argument('-m2', '--matrix2', type=str, required=True, help="gene expression per cell for species 1")

PARSER.add_argument('-c1', '--clusters1', type=str, required=True, help="cell-to-clusters table for species 1")

PARSER.add_argument('-c2', '--clusters2', type=str, required=True, help="cell-to-clusters table for species 2")

PARSER.add_argument('-f', '--gene_families', type=str, required=True, help="gene families")

PARSER.add_argument('-c', '--cores', type=int, required=False, default=1)

PARSER.add_argument('-o', '--output_dir', type=str, required=False, default="out_scCompare/")

PARSER.add_argument('-sp1', '--label_species1', type=str, required=False, default="")

PARSER.add_argument('-sp2', '--label_species2', type=str, required=False, default="")

PARSER.add_argument('--filter_out', type=str, required=False, default="")

PARSER.add_argument('--random', action='store_true', help="this triggers a randomization of the gene homologies to evaluate the expected random background correlations")

# PARSER.add_argument('--reoptimize', action='store_true', help="this triggers a re-optimization of ec scores on the whole set after paralogs selection")

PARSER.add_argument('--no_paralogs', action='store_true', help="only use 1-to-1 orthologs")

ARGS = vars(PARSER.parse_args())

# Create output dir
ARGS['output_dir'] = ARGS['output_dir'].strip('/') + '/' 
os.makedirs(ARGS['output_dir']+ 'random', exist_ok=True)

print('# Matrix 1 #')
norm_mat1 = ARGS['output_dir'] + Path(ARGS['matrix1']).stem + '_expr_clusters_norm.tsv'
nm.normalize_main(ARGS['matrix1'], ARGS['clusters1'], norm_mat1, cores=ARGS['cores'], filter_out_start=ARGS['filter_out'])

print('# Matrix 2 #')
norm_mat2 = ARGS['output_dir'] + Path(ARGS['matrix2']).stem + '_expr_clusters_norm.tsv'
nm.normalize_main(ARGS['matrix2'], ARGS['clusters2'], norm_mat2, filter_out_start=ARGS['filter_out'])

print('# Expression conservation of 1-1 orthologs + paralogs selection #')
sp.select_paralogs_main(norm_mat1, norm_mat2, ARGS['gene_families'], ARGS['output_dir'], max_combin=200,
                        ncores=ARGS['cores'], noparalogs=ARGS['no_paralogs'])

print('# Correlation between clusters of species 1 and clusters of species 2 #')
cp.compare_main(norm_mat1, norm_mat2, ARGS['output_dir'], ARGS['cores'], ARGS['label_species1'], ARGS['label_species2']) #, random=random_homologs_dir reoptimize_ec=ARGS['reoptimize'] #remove the reoptimize arg

print('Done!')

"""
python -m scCompare --matrix1 ../data/Cg_matrix_EM.csv --matrix2 ../data/Pc_matrix_EM.csv --clusters1 ../data/Cragig_cell_id.tsv --clusters2 ../data/Procro_cell_id.tsv --gene_families ../data/orthologous_pairs.txt -c 20 -sp1 Cragig -sp2 Procro
"""

"""
python -m scCompare --matrix1 ../data/placozoa/Tadh_final_matrix.tsv --matrix2 ../data/placozoa/TrH2_final_matrix.tsv --clusters1 ../data/placozoa/Tadh.sc_annot.tsv --clusters2 ../data/placozoa/TrH2.sc_annot.tsv --gene_families ../data/placozoa/orthologous_pairs_ok.txt -c 20 -sp1 Tadh -sp2 TrH2
"""