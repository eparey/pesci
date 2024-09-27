# pesci - Pretty Easy Single-cell Comparisons using ICC

[![python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm applied to the **comparison of single-cell datasets across two study species** ([Najle, Grau-Bové et al.](https://doi.org/10.1016%2Fj.cell.2023.08.027)).

Pesci takes as input **single-cell expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies** files.


## Citing
	
- Najle, Grau-Bové et al. Stepwise emergence of the neuronal gene expression program in early animal evolution. (2023). Cell.

- pesci application note


## Installation

- With conda (local for now, will upload to bioconda at a later stage)

	`git clone git@github.com:eparey/pesci.git`

	`cd pesci`

	`mamba install pesci -c ./recipe/build -c conda-forge`


- With pip (local for now, will upload to pypi at a later stage)

	`git clone git@github.com:eparey/pesci.git`

	`cd pesci`

	`pip install . -I` #local install with pip, suggest installing in an isolated env (conda or venv)


## Usage

For a description of accepted input formats and available options, please refer to `pesci --help` or read the full documentation [TODO link].

- To run pesci on provided example data (datasets from Piovani et al., 2023) [TODO link]:

```
pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -sp1 Oyster-larva -sp2 Flatworm-larva --colbroad broad --cores 4
```

 ![pesci fig](https://github.com/eparey/pesci/blob/main/docs/img/Oyster-larva-Flatworm-larva_correlation_scores_matrix.png)


## License

## Contacts