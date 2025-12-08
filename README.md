# pesci - Pretty Easy Single-cell Comparisons using ICC

![Static Badge](https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9%7C3.10%7C3.11%7C3.12-blue?logo=python)

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm applied to the **comparison of single-cell gene expression datasets across two study species**, first proposed by ([Najle, Grau-Bové et al.](https://doi.org/10.1016%2Fj.cell.2023.08.027)).


## Citing

- pesci application note

- Najle, Grau-Bové et al. Stepwise emergence of the neuronal gene expression program in early animal evolution. (2023). Cell.


## Installation

- With conda (local for now, will upload to bioconda at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`conda create -n pesci_env python=3.12 && conda activate pesci_env` #suggest installing in an isolated env 'pesci_env' (conda or venv)

	`conda install pesci -c ./recipes/build -c conda-forge`


- With pip (local for now, will upload to pypi at a later stage)

	`git clone git@github.com:eparey/pesci.git && cd pesci`

	`conda create -n pesci_env python=3.12 && conda activate pesci_env` #suggest installing in an isolated env 'pesci_env' (conda or venv)

	`pip install . -I`


## Usage

Pesci takes as input **single-cell expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies** files.

```
pesci -m1 mat_sp1.h5ad -m2 mat_sp2.tsv -c1 cell_id_sp1.tsv -c2 cell_id_sp2.tsv -g orthologs.txt 
```

For a description of accepted input formats and available options, please refer to `pesci --help`.

- To run pesci on provided example data (datasets from [Piovani et al., 2023](https://doi.org/10.1126/sciadv.adg6034)):

```
pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Procro_matrix_EM.tsv.gz -c1 data/Cragig_matrix_EM.tsv.gz -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -sp1 Oyster-larva -sp2 Flatworm-larva --colbroad broad --cores 4
```

![pesci fig](https://github.com/eparey/pesci/blob/main/docs/img/Oyster-larva-Flatworm-larva_correlation_scores_matrix.png)


## Contacts

- [Elise Parey](e.parey@ucl.ac.uk)
- [Laura Piovani](l.piovani@ucl.ac.uk)

