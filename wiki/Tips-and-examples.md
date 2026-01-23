# Tips and Example Gallery

This page presents useful tips and illustrated examples of different potential use cases of Pesci.

## Table of Contents
1. [Customizing the output figure](#customizing-the-output-figure)
2. [Providing multiple libraries as input](#providing-multiple-libraries-as-input)
3. [Exploring the impact of the random seed](#exploring-the-impact-of-the-random-seed)
4. [Using a subset of the cell clusters](#using-a-subset-of-the-cell-clusters)
5. [Comparing organs within a species](#comparing-organs-within-a-species)


## Customizing the output figure

Different options are available to customize the output heatmap figure. For operations not covered by available options, Pesci also saves the underlying matrix to a file that can be loaded into python and R.

### Default

The default is to attempt to maximize 1-1 cell cluster matches on the diagonal, this works well for closely related species, but is less suitable to our example use case:

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_default.png)

### Reordering the rows and columns with --reorder

The option `--reorder` accepts the following arguments to reorder rows and columns: `Diag` (default, see above), `Clust` uses hierarchical clustering (average linkage, correlation), and `Alpha` for alphabetical: 

#### Clust

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder Clust
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_clust.png)

#### Alpha

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder None
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_none.png)


### Ordering by broad annotation

Alternatively, if broad annotation are supplied in addition to the cell cluster labels (as an additional column in the cell-to-cluster file or the .h5ad scanpy object), these will be used for columns and rows reordering in the output heatmap:

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad.png)


### Adding species labels

Providing species labels to -l1 and -l2 is not only useful to trace the output files of a specific pesci run, it also labels the heatmap axes:

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -l1 Oyster-larva -l2 Flatworm-larva
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad_wnames.png)


### Color palette

The continuous color palette used for the heatmap can be changed useing the `--seaborn_cmap` argument, which accepts any named [Color Brewer](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html) palette.

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -l1 Oyster-larva -l2 Flatworm-larva --seaborn_cmap Blues
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad_wnames_cmap.png)

### Advanced 

Finally, the raw results matrix (`output_pesci/Oyster-larva-Flatworm-larva_correlation_scores_matrix.csv`) can be directly loaded into python or R for more direct customization.

## Providing multiple libraries as input

Pesci supports use cases where single-cell expression datasets were generated from samples and where these haven't been merged into a single file.

There is one important requirement for these files to be successfully loaded by pesci: cell barcodes should be distinct across input data matrices (consider adding a sample prefix to cell barcodes if necessary), as they correspond to distinct cells. For cases where the same library was sequenced multiple times, these should be merged into a single file prior to running pesci. 

For inputs requiring a cell-to-cluster correspondence file, this file should remain as a single file containing all of the cell barcodes.

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

- To turn-off down-sampling (slower):
	```
	pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --do_not_downsample  --cores 4 --force
	```

Pesci reports the total number of one-to-one orthologs included in the comparisons (e.g. `[INFO]: Found 4567 one-to-one orthologs`) and total number of orthologs considered after inclusion of many-to-many and many-to-one (e.g. `Total retained genes for cross-species comparison: 6222 genes.`). The number of total included genes has an important impact on the cell matches found: the more genes are included, the more likely it is that cell markers are included, hence increasing expression similarity scores. It is therefore important to provide comprehensive orthologies to pesci.


## Using a subset of the cell clusters

Two options (`--filter_out` and `--keep_only`) are available to restrict comparison to only a subset of the cell clusters, for instance to exclude artefactual or poorly characterized cell clusters, or to compare only specific cell types.

These options take one or several string of characters, comma-separated, excluding any clusters whose name starts with the provided strings. While `--filter_out` and `--keep_only` apply to both species, `--filter_out1` and `--keep_only1` (respectively `--filter_out2` and `--keep_only2` ) apply to species 1 only (resp. species 2).


- Example 1: filtering out all clusters starting with "Unk" in both species:

	```
	pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --filter_out 'Unk' --colbroad broad -l1 Oyster-larva -l2 Flatworm-larva
	```
![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrixUnk.png)


- Example 2: filtering out all clusters starting with "Unk,1,6,8,9" species2:


	```
	pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --filter_out2 'Unk,1,6,8,9' --colbroad broad -l1 Oyster-larva -l2 Flatworm-larva
	```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrixfilttflat.png)

- Example 3: keeping only bipolar

example (for instance ferret - cow) bp 

example (for instance ferret - cow) bp + marker specificity




## Comparing organs within a species

Pesci can be run to compare two single-cell expression datasets from different organs but of the same species, by using the `--within_species` file. Here, the orthology file is obviously not necessary and can be omitted, identical gene names will be automatically matched across the input matrices.

- Example: comparison of fruit fly trachea and proboscis from the [FLY CELL ATLAS](https://flycellatlas.org/):

```
pesci --matrix1 FlyCellAtlas_adult_trachea_sparse.h5ad  --matrix2 FlyCellAtlas_adult_proboscis_sparse.h5ad  --clusters1 'annotation' --clusters2 'annotation' -o pesci_fly --within_species --filter_out 'male,ovary' --seaborn_cmap Purples -l1 "Adult-Fly-Trachea" -l2 "Adult-Fly-Proboscis"
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_fly.png)

