# pesci - Pretty Easy Single-cell Comparisons with ICC

TODO badges

TODO LOGO

TODO short description

TODO HOWTOCITE (placozoa paper + pesci application note)

## Quick start

### Installation


#### Conda (local for now, will upload to bioconda at a later stage)

`git clone git@github.com:eparey/pesci.git`

`cd pesci`

`mamba install pesci -c ./recipe/build -c conda-forge`


#### Pip (local for now, will upload to pypi at a later stage)

`git clone git@github.com:eparey/pesci.git`

`cd pesci`

`pip install . -I` #local install with pip, suggest installing in an isolated env (conda or venv)


### Usage

Test the installation on provided example data:

`pesci -m1 data/Cg_matrix_EM.tsv.gz -m2 data/Pc_matrix_EM.tsv.gz -c1 data/Cragig_cell_id.tsv -c2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt --cores 4 -sp1 Cragig -sp2 Procro`

TODO Link to the full doc for info on all options (see --help)

TODO License