# Tips and Example Gallery

This page presents useful tips and illustrated examples of different potential use cases of Pesci.

## Table of Contents
1. [Providing multiple libraries as input](#providing-multiple-libraries-as-input)
2. [Exploring the impact of the random seed](#exploring-the-impact-of-the-random-seed)
3. [Using a subset of the cell clusters](#using-a-subset-of-the-cell-clusters)
4. [Comparing organs within a species](#comparing-organs-within-a-species)
5. [Customizing the output figure](#customizing-the-output-figure)


## Providing multiple libraries as input

Pesci supports use cases where single-cell expression datasets were generated from multiple sequenced libraries and where these haven't been merged into a single file.

There is one important requirement for these files to be successfully loaded by pesci: cell barcodes should be distinct across input data matrices (consider adding a sample prefix to cell barcodes if necessary). Since these correspond to distinct single-cell RNA-seq libraries, they contain distinct cells. For cases where the same library was sequenced multiple times, these should be merged into a single file prior to running pesci. 

Note also that for inputs requiring a cell-to-cluster correspondence file, this should be a single file containing all of the cell barcodes.

- Provided these conditions are met, several input matrix files, in any of the accepted formats, can be supplied to `--matrix1` (or `--matrix2`), separated by spaces:

	```
	pesci --matrix1 data/Cg_matrix_EM_part1.tsv.gz data/Cg_matrix_EM_part2.tsv.gz data/Cg_matrix_EM_part3.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --cores 4
	```


## Exploring the impact of the random seed

Two options are available to explore the impact of the random seed in the observed results: `--seed` and `--do_not_downsample`. For computational efficiency, Pesci randomly samples 1000 one-to-one orthologous genes to be used as a "skeleton" to select best pairs for the many-to-many / many-to-one ortholog cases. In our hands, this had little impact on the results, but is relevant to explore on different datasets.

- To change the random seed:

	```
	pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --seed 666  --cores 4 --force
	```

- To turn-off down-sampling:
	```
	pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --do_not_downsample  --cores 4 --force
	```

Pesci reports the total number of one-to-one orthologs included in the comparisons (e.g. `[INFO]: Found 4567 one-to-one orthologs`) and total number of orthologs considered after inclusion of many-to-many and many-to-one (e.g. `Total retained genes for cross-species comparison: 6222 genes.`). The number of total included genes has an important impact on the cell matches found: the more genes are included, the more likely it is that cell markers are included, hence increasing expression similarity scores. It is therefore important to provide comprehensive orthologies to pesci.


## Using a subset of the cell clusters

Two options (`--filter_out` and `--keep_only`) are available to restrict comparison to only a subset of the cell clusters, for instance to exclude artefactual or poorly characterized cell clusters, or to compare only specific cell types.

example filter out unknown clusters in sp1 

example filter out unknown clusters in both


example (for instance ferret - cow) bp 

example (for instance ferret - cow) bp + marker specificity




## Comparing organs within a species

same species but different organs? --> illustrated example of the droso 



## Customizing the output figure

Different options are available to customize the output heatmap figure. For operations not covered by available options, Pesci also saves the underlying matrix to a file that can be loaded into python and R.

default

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_default.png)


--reorder Clust
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder Clust
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_clust.png)


--reorder None
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder None
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_none.png)


using broad annotation 
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad.png)

add species labels
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -sp1 Oyster-larva -sp2 Flatworm-larva
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad_wnames.png)

changing the palette
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -sp1 Oyster-larva -sp2 Flatworm-larva --seaborn_cmap Blues
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad_wnames_cmap.png)

Finally, the raw results matrix (`output_pesci/Oyster-larva-Flatworm-larva_correlation_scores_matrix.csv`) can be directly loaded into python or R for more direct customization.
