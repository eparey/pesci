# pesci - Pretty Easy Single-cell Comparisons with ICC

TODO badges

TODO LOGO

TODO short description

TODO HOWTOCITE

## Quick start

### Installation

- Install pesci

`git clone git@github.com:eparey/pesci.git`
`cd pesci`
`pip install . -I` #local install with pip, suggest installing in an isolated env (conda or venv)

- To get the example data and run a test job:

`git lfs install` (if you do not have lfs, see instructions [here](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage))
`git lfs pull`

### Usage

`pesci -m1 data/Cg_matrix_EM.tsv -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv.gz -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt -c 8`

TODO Link to the full doc for info on all options (see --help)

TODO License