# pesci - Pretty Easy Single-cell Comparisons using ICC

[![python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm applied to the **comparison of single-cell datasets across two study species** ([Najle, Grau-Bové et al.](https://doi.org/10.1016%2Fj.cell.2023.08.027)).

Pesci takes as input **single-cell expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies** files.


## Citing
	- Najle, Grau-Bové et al. (2023). Stepwise emergence of the neuronal gene expression program in early animal evolution. Cell 

	- pesci application note


## Installation

#### Conda (local for now, will upload to bioconda at a later stage)

`git clone git@github.com:eparey/pesci.git`

`cd pesci`

`mamba install pesci -c ./recipe/build -c conda-forge`


#### Pip (local for now, will upload to pypi at a later stage)

`git clone git@github.com:eparey/pesci.git`

`cd pesci`

`pip install . -I` #local install with pip, suggest installing in an isolated env (conda or venv)


## Usage

Test the installation on provided example data:

```
pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -c 4 -sp1 Oyster-larva -sp2 Flatworm-larva
```

TODO Figure

Use broad annotation to order the heatmap:
```
pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -c 4 -sp1 Oyster-larva -sp2 Flatworm-larva --force --colbroad 
```

TODO Figure

TODO Link to the full doc for info on all options (see --help)

TODO License