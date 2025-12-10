List of available options
==========================

This page provides a description of all available options when running pesci. For more details on input and output files, including R code to prepare input files from Seurat objects, please see `Help With Input Files <https://github.com/eparey/pesci/blob/main/wiki/Input-Files.md>`_ and `Understanding Output files <https://github.com/eparey/pesci/blob/main/wiki/Outputss.md>`_. For example use cases with specific options, please see `Pesci Examples <https://github.com/eparey/pesci/blob/main/wiki/Examples.md>`_.


Minimal command-line examples
-----------------------------

Below are minimal examples to run pesci, each showcasing different accepted formats for input single-cell expression datasets. The corresponding data files are available as examples in the `data/` directory.

- **Example 1**: species 1 data as a sparse count matrix (same format as a cellranger directory) and species 2 provided as a dense count matrix (least optimal format)

.. code-block:: sh

  pesci --matrix1 data/Cragig_sparse_data/  --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt

- **Example 2**: species 1 data as a scanpy h5ad and species 2 provided as a dense count matrix (least optimal format) - note that here 'cluster_name' is the name of the column with cluster annotation in the h5ad

.. code-block:: sh

  pesci --matrix1 data/Cragig_matrix.h5ad  --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 'cluster_name' --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt

- **Example 3**: data for both species provided as dense count matrices, but with species 1 data from 3 different libraries with counts respectively saved in 3 different files

.. code-block:: sh

  pesci --matrix1 data/Cg_matrix_EM_part1.tsv data/Cg_matrix_EM_part2.tsv data/Cg_matrix_EM_part3.tsv --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt



note 

files can be compressed in .gz or not.

*important* note somewhere: if running different comparisons use sp. labels and/or different output directory, otherwise with warning will reuse saying already computed


Usage
------

All available options are described below and can be printed using `pesci --help`. Note that flags (--force, --ono2one_only, --do_not_plot_warn, --show_auto_threshold, --no_pbar) do not require any arguments to be set, for instance adding`--force` will effectively force recomputation of all intermediary files.

**Command-Line Options (Markdown Table Format)**



## **Options**
| Flags | Argument | Description |
|-------|----------|-------------|
| `-h`, `--help` | — | Show help message and exit. |
| `-v`, `--version` | — | Show program version and exit. |

---

## **Required Arguments**
| Flags | Argument | Description |
|-------|----------|-------------|
| `--matrix1` | `MATRIX1` | Gene expression per cell for species 1 (raw count matrix). Accepts: dense text (optionally `.gz`), Cell Ranger directory, or `.h5ad`. Multiple inputs allowed. |
| `--matrix2` | `MATRIX2` | Gene expression per cell for species 2. Same format options as `--matrix1`. Multiple inputs allowed. |
| `--clusters1` | `CLUSTERS1` | Cell-to-clusters table for species 1, or cluster column name in `.h5ad`. |
| `--clusters2` | `CLUSTERS2` | Cell-to-clusters table for species 2, or cluster column name in `.h5ad`. |
| `-g,⁠--orthologous_gene_pairs` | `ORTHOLOGOUS_GENE_PAIRS` | Tab-delimited orthology file (one gene pair per line). |
---

## **Recommended Arguments**
| Flags | Argument | Description |
|-------|----------|-------------|
| `-o`, `--outdir` | `OUTDIR` | Output directory (created if missing). **Default:** `output_pesci/` |
| `--label_species1` | `LABEL_SPECIES1` | Species 1 name used in plots/outputs. **Default:** `sp1` |
| `--label_species2` | `LABEL_SPECIES2` | Species 2 name used in plots/outputs. **Default:** `sp2` |

---

## **Resources**
| Flags | Argument | Description |
|-------|----------|-------------|
| `-c`, `--cores` | `CORES` | Number of CPU cores for parallel operations. **Default:** `1` |

---

## **Running**
| Flags | Argument | Description |
|-------|----------|-------------|
| `--force` | — | Recompute all intermediates even if they exist. **Default:** `False` |
| `--seed` | `SEED` | Random seed. **Default:** `123` |
| `--marker_specificity` | `MARKER_SPECIFICITY` | Marker specificity (0.5–0.95). Controls how fold change is computed (median vs percentile). **Default:** `0.5` |
| `--ono2one_only` | — | Use only 1-to-1 orthologs. **Default:** `False` |
| `--ec_threshold_many` | `EC_THRESHOLD_MANY` | Retain all distinct high-score pairs (> threshold) for many-to-many and many-to-one orthologs. **Default:** `None` |
| `--do_not_downsample` | — | Use all 1-to-1 orthologs (no random subsampling). **Default:** `False` |
| `--random_id` | `RANDOM_ID` | Randomizes orthologies and stores results under `outdir/random_id`. **Default:** blank |
| `--no_pbar` | — | Disable progress bar. **Default:** `False` |

---

## **Inputs**
| Flags | Argument | Description |
|-------|----------|-------------|
| `--min_umi` | `MIN_UMI` | Minimum total UMI per gene to retain. **Default:** `10` |
| `--colclust` | `COLCLUST` | Column containing cluster labels (ignored for `.h5ad`). Defaults to `cluster_name`, else second column. |
| `--colclust1` | `COLCLUST1` | Same as `--colclust` but for species 1. **Default:** `None` |
| `--colclust2` | `COLCLUST2` | Same as `--colclust` but for species 2. **Default:** `None` |
| `--colbroad` | `COLBROAD` | Column for broad annotations; used only for heatmap. **Default:** `None` |
| `--colbroad1` | `COLBROAD1` | Broad annotations column for species 1. **Default:** `None` |
| `--colbroad2` | `COLBROAD2` | Broad annotations column for species 2. **Default:** `None` |
| `--filter_out` | `FILTER_OUT` | Remove clusters starting with given prefixes (comma-separated). Applies to both species. |
| `--filter_out1` | `FILTER_OUT1` | Same as above, species 1 only. |
| `--filter_out2` | `FILTER_OUT2` | Same as above, species 2 only. |
| `--keep_only` | `KEEP_ONLY` | Keep only clusters matching given prefixes (comma-separated). Applies to both species. |
| `--keep_only1` | `KEEP_ONLY1` | Same as above, species 1 only. |
| `--keep_only2` | `KEEP_ONLY2` | Same as above, species 2 only. |

---

## **Outputs**
| Flags | Argument | Description |
|-------|----------|-------------|
| `-l`, `--logfile` | `LOGFILE` | Write log to this file instead of default (`pesci.log`). |
| `-f`, `--figure_format` | `FIGURE_FORMAT` | Output figure format: `svg`, `png`, `pdf`. **Default:** `pdf` |
| `-r`, `--reorder` | `REORDER` | Ordering method for heatmap clusters: `DiagKeep`, `Clust`, or `None`. Ignored with broad annotations. **Default:** `DiagKeep` |
| `--min_fc` | `MIN_FC` | Minimum fold change to count a marker gene (used for co-expressed marker table). **Default:** `1.5` |
| `--seaborn_cmap` | `SEABORN_CMAP` | Seaborn colormap for heatmaps. **Default:** `BuPu` |
| `--show_auto_threshold` | — | Highlight matches above thresholds on heatmap. **Default:** `False` |
| `--do_not_plot_warn` | — | Do not display heatmap warning for matches supported by <10 genes. **Default:** `True` |
