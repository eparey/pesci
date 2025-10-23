# pesci - Pretty Easy Single-cell Comparisons using ICC

![Static Badge](https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9%7C3.10%7C3.11%7C3.12-blue?logo=python)

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm applied to the **comparison of single-cell gene expression datasets across two study species**, first proposed by ([Najle, Grau-Bové et al.](https://doi.org/10.1016%2Fj.cell.2023.08.027)).


## Citing

- pesci application note

- Najle, Grau-Bové et al. Stepwise emergence of the neuronal gene expression program in early animal evolution. (2023). Cell.


## Installation

- With conda (local for now, will upload to bioconda at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`conda install pesci -c ./recipes/build -c conda-forge`


- With pip (local for now, will upload to pypi at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`pip install . -I` #local install with pip, suggest installing in an isolated env (conda or venv)


## Usage

Pesci takes as input **single-cell expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies** files.

```
pesci -m1 mat_sp1.h5ad -m2 mat_sp2.tsv -c1 cell_id_sp1.tsv -c2 cell_id_sp2.tsv -g orthologs.txt 
```

For a description of accepted input formats and available options, please refer to `pesci --help`.

- To run pesci on provided example data (datasets from [Piovani et al., 2023](https://doi.org/10.1126/sciadv.adg6034)):

```
pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -sp1 Oyster-larva -sp2 Flatworm-larva --colbroad broad --cores 4
```

![pesci fig](https://github.com/eparey/pesci/blob/main/docs/img/Oyster-larva-Flatworm-larva_correlation_scores_matrix.png)


## License

## Contacts

## TODO: make a more complete doc with Gihtub Pages and link to it

For instance, one way to format from seurat:

```
library(Matrix)
library(R.utils)
library(data.table)
library(tidyverse)

counts <- mySeuratObj@assays$RNA@counts
# Output counts matrix
writeMM(counts, paste0(data_dir, 'matrix.mtx'))
gzip(paste0(data_dir, 'matrix.mtx'))

# Output cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes), paste0(data_dir, 'barcodes.tsv'),
           col_names = FALSE)
gzip(paste0(data_dir, 'barcodes.tsv'))

# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", paste0(data_dir, 'features.tsv'),
           col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))
```

NOTES TODO: important in data preperation to have same gene id in the orthology and matrix files AND that gene ids are unique for each species (recommend adding prefix with sp name)