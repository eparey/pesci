# pesci - Pretty Easy Single-cell Comparisons using ICC

![Static Badge](https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9%7C3.10%7C3.11%7C3.12-blue?logo=python)

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm to compare **single-cell gene expression datasets across or within species**. 

The ICC algorithm was first proposed by [Tirosh and Barkai](https://doi.org/10.1186/gb-2007-8-4-r50) and first applied to single-cell gene expression by ([Najle, Grau-Bové et al.](https://doi.org/10.1016%2Fj.cell.2023.08.027)).


## Citing

- Pesci -- Coming soon

- Najle, Grau-Bové et al. Stepwise emergence of the neuronal gene expression program in early animal evolution. (2023). Cell.


## Installation

- With conda (local for now, will upload to bioconda at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`conda create -n pesci_env python=3.12 && conda activate pesci_env` #suggest installing in an isolated env 'pesci_env' (conda or venv) 

	`conda install pesci -c ./recipes/build -c conda-forge`


- With pip (local for now, will upload to pypi at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`conda create -n pesci_env python=3.12 && conda activate pesci_env` #suggest installing in an isolated env 'pesci_env' (conda or venv) [may need python 3.9 to install llvmlite dependency on macos...]

	`pip install . -I`


## Usage

Pesci takes as input **single-cell gene expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies**.

```sh
pesci --matrix1 mat_sp1.h5ad --matrix2 mat_sp2.tsv --clusters1 cluster_annotations --clusters2 cell_id_sp2.tsv --ortho_pairs orthologs.txt 
```

For a description of accepted input formats and available options, please refer to `pesci --help`.

- To run pesci on provided example data (datasets from [Piovani et al., 2023](https://doi.org/10.1126/sciadv.adg6034)):

```sh
pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv --ortho_pairs data/orthologous_pairs_Procro-Cragig.txt -l1 Oyster-larva -l2 Flatworm-larva --colbroad broad --cores 4
```

![pesci fig](https://github.com/eparey/pesci/blob/main/wiki/img/Oyster-larva-Flatworm-larva_correlation_scores_matrix.png)


## Documentation

LINK TO WIKI HERE TODO

## Contacts

- [Elise Parey](e.parey@ucl.ac.uk)
- [Laura Piovani](l.piovani@ucl.ac.uk)

