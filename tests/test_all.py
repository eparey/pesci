"""
Tests for pesci cli.

This should be improved with unit tests, this only checks that the code does not get broken from 
the first working version (0.1.0) onwards.

Only default parameters and input as raw count matrix are covered.

"""

import os
import pytest
from cli_test_helpers import shell

import filecmp

import pesci

import pesci.normalize as nm
import pesci.iterative_comparison_coexpression as icc
import pesci.compare as cp


def test_install():
    """
    check install
    """
    result = shell("pesci --help")
    assert result.exit_code == 0


@pytest.mark.long
def test_normalize():
	"""
	check normalize (default parameters) against the oyster data
	"""
	mat1 = 'data/Cg_matrix_EM.tsv.gz'
	cl1 = 'data/Cragig_cell_id.tsv'
	output = 'tests/normed_tmp.tsv'
	expected_out = 'data/tests/Cragig_Cg_matrix_EM_expr_clusters_norm_all.tsv'

	#test normalize function with Cragig
	nm.normalize([mat1], cl1, output)

	#check results against expected
	assert filecmp.cmp(output, expected_out)
	os.remove(output)


@pytest.mark.long
def test_icc():
	"""
	check icc implementation against the oyster-flatworm data
	"""
	norm_mat1 = 'data/tests/Cragig_Cg_matrix_EM_expr_clusters_norm_all.tsv'
	norm_mat2 = 'data/tests/Procro_Pc_matrix_EM_expr_clusters_norm_all.tsv'
	ortho = 'data/orthologous_pairs_Procro-Cragig.txt'

	#test icc by loading Cragig - Procro normed
	icc.icc(norm_mat1, norm_mat2, ortho, 'tests/sp1-sp2', ncores=4, seed=123)

	expected_out = 'data/tests/files/Cragig-Procro_1-to-1-orthologs_correlation_scores.tsv'
	output = 'tests/sp1-sp2_1-to-1-orthologs_correlation_scores.tsv'
	assert filecmp.cmp(output, expected_out)
	os.remove(output)

	output = 'tests/sp1-sp2_orthologs_many_correlation_scores.txt'
	expected_out = 'data/tests/files/Cragig-Procro_orthologs_many_correlation_scores.txt'
	assert filecmp.cmp(output, expected_out)
	os.remove(output)

	os.remove('tests/sp1-sp2_skipped_ogs.txt')


@pytest.mark.long
def test_compare():
	"""
	check cell clusters comparisons against the oyster-flatworm data
	"""
	#we have an issue for tests here: input will go in same dir as output
	norm_mat1 = 'data/tests/Cragig_Cg_matrix_EM_expr_clusters_norm_all.tsv'
	norm_mat2 = 'data/tests/Procro_Pc_matrix_EM_expr_clusters_norm_all.tsv'
	

	#test compare by loading icc results and Cragig - Procro normed and ec scores
	cp.compare(norm_mat1, norm_mat2, 'data/tests/', sp1='Cragig', sp2='Procro')
	output = 'data/tests/Cragig-Procro_correlation_scores_matrix.csv'
	expected_out = 'data/tests/result_correlation_scores_matrix.csv'
	assert filecmp.cmp(output, expected_out)
	os.remove(output)

	os.remove('data/tests/Cragig-Procro_correlation_scores_matrix.svg')
	os.remove('data/tests/Cragig-Procro_expression_conservation_scores.csv')
	os.remove('data/tests/Cragig-Procro_gene_coexpression_table.csv')

	os.remove('data/tests/files/Cragig-Procro_best_match.txt')
