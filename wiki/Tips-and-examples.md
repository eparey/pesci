# Tips and Example Gallery

This page presents useful tips and illustrated examples of different potential use cases of Pesci.

## Several lib.
several libraries were sequenced and stored in distinct files


## Impact of the random seed

impact of random seed --> show she-mes match goes away depending on seed (can also use all orthologs)
paying attention to the number of orthologs, importance of marker specificity


## Customizing the output figure

Different options are available to customize the output heatmap figure. For operations not covered by available options, Pesci also saves the underlying matrix to a file that can be loaded into python and R.

default

```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt
```

--reorder Clust
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder Clust
```

--reorder None
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --reorder None
```


using broad annotation 
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad
```

add species labels
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -sp1 Oyster-larva -sp2 Flatworm-larva
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/matrix_broad_wnames.png)

changing the palette
```
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt --colbroad broad -sp1 Oyster-larva -sp2 Flatworm-larva ----seaborn_cmap Blues
```

(can also work on the raw result .tsv)

Loading in python
python```
test
```

Loading in R
R```
test
```

## Sub-setting

example closer distance comparison + sub-setting neuronal cell types (for instance ferret - cow)

--filter_out and --keep_only


## Using pesci to compare datasets from the same species

same species but different organs? --> illustrated example of the droso 

