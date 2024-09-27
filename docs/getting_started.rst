Quick start
=============

Installation
-------------
.. prompt:: bash

	conda install -c conda-forge -c bioconda pesci

or

.. prompt:: bash

	pip install pesci

Command-line example
---------------------

.. prompt:: bash

	pesci --matrix1 mat1.tsv.gz --matrix2 mat2_rep1.tsv mat2_rep2.tsv -clusters1 cell_id1.tsv --clusters2 cell_id2.tsv -g orthologous_pairs.txt


Usage
------

**options:**
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

**required arguments:**
  --matrix1 MATRIX1
                        gene expression per cell for species 1 (raw count matrix), either: a dense matrix text file (can be .gz), a cellranger directory or an h5ad file; several inputs can be specified

  --matrix2 MATRIX2
  						gene expression per cell for species 2 (raw count matrix), either: a dense matrix text file (can be .gz), a cellranger directory or an h5ad file; several inputs can be specified

  --clusters1 CLUSTERS1
                        cell-to-clusters table for species 1 or name of cluster column in h5ad

  --clusters2 CLUSTERS2
                        cell-to-clusters table for species 2 or name of cluster column in h5ad

  -g ORTHOLOGOUS_GENE_PAIRS, --orthologous_gene_pairs ORTHOLOGOUS_GENE_PAIRS
                        gene orthologies file (tab-delimited, 1 pair per line)

**resources:**
  -c CORES, --cores CORES
                        number of cores to use for parallel operations (default: 1)

**running:**
  --force               recompute the per cluster normalized gene expression and ec scores even if output already exists (default: False)
  --seed SEED           set random seed (default: 123)
  --marker_specificity MARKER_SPECIFICITY
                        marker specificity (btwn 0.5 and 0.95) - used in fold change expression calculation
                        as follows: 0.5 computes fold change in cluster 'a' vs median expression across
                        clusters, 0.75 computes fold change in cluster 'a' vs 75 expression percentile across
                        clusters (3rd quartile) (default: 0.5)
  --ono2one_only        only use 1-to-1 orthologs (default: False)
  --ec_threshold_many EC_THRESHOLD_MANY
                        if set, instead of only best pair retain all distinct pairs with score > threshold
                        (default: None)
  --random_id RANDOM_ID
                        triggers a randomization of orthologies and store results in outdir/random_id
                        (default: )

**input:**
  --min_umi MIN_UMI     Minimum total umi for a gene to be retained for comparison (default: 10)
  --colclust COLCLUST   name of the column corresponding to clusters in the cell to cluster annotation file
                        (does not apply to h5ad input as it is specified by --clusters1/--cluster2 directly)
                        - if not set, will use column `cluster_name`, if it does not exist, will attempt to
                        use the 2nd column

                        use --colclust1 & --colclust2 to specify different column names
                        for species1 & species2. (default: cluster_name)
  --colclust1 COLCLUST1
                        name of the column corresponding to clusters in the cell to cluster annotation file
                        for species1 (does not apply to h5ad input as it is specified by --clusters1
                        directly) - if not set, will use column `cluster_name`, if it does not exist, will
                        attempt to use the 2nd column. (default: None)
  --colclust2 COLCLUST2
                        name of the column corresponding to clusters in the cell to cluster annotation file
                        for species2 (does not apply to h5ad input as it is specified by --clusters1
                        directly) - if not set, will use column `cluster_name`, if it does not exist, will
                        attempt to use the 2nd column. (default: None)
  --colbroad COLBROAD   name of the column corresponding to broad annotations in the cell to cluster
                        annotation file or h5ad - if not set, will not use broad annotations (this is for the
                        heatmap plot only)

                        use --colbroad1 & --colbroad2 to specify different column names
                        for species1 & species2. (default: None)
  --colbroad1 COLBROAD1
                        name of the column corresponding to broad annotations in the cell to cluster
                        annotation file or h5ad of species1 - if not set, will not use broad annotations
                        (this is for the heatmap plot only). (default: None)
  --colbroad2 COLBROAD2
                        name of the column corresponding to broad annotations in the cell to cluster
                        annotation file or h5ad of species2 - if not set, will not use broad annotations
                        (this is for the heatmap plot only). (default: None)
  --filter_out FILTER_OUT
                        discard clusters starting with specified string,several can be specified as a comma-
                        separated list (ex: --filter_out BC\_,RGC\_). Applies to both species1 & species2.
                        (default: )
  --filter_out1 FILTER_OUT1
                        discard clusters starting with specified string,several can be specified as a comma-
                        separated list (ex: --filter_out1 BC\_,RGC\_) - applies to species1 only. (default: )
  --filter_out2 FILTER_OUT2
                        discard clusters starting with specified string,several can be specified as a comma-
                        separated list (ex: --filter_out2 BC\_,RGC\_) - applies to species2 only. (default: )
  --keep_only KEEP_ONLY
                        keep only clusters starting with specified string, several can be specified as a
                        comma-separated list (ex: --keep_only BC\_,RGC\_) - applies to both species1 &
                        species2. (default: )
  --keep_only1 KEEP_ONLY1
                        keep only clusters starting with specified string, several can be specified as a
                        comma-separated list (ex: --keep_only1 BC\_,RGC\_) - applies to species1 only.
                        (default: )
  --keep_only2 KEEP_ONLY2
                        keep only clusters starting with specified string, several can be specified as a
                        comma-separated list (ex: --keep_only2 BC\_,RGC\_) - applies to species2 only.
                        (default: )

**output:**
  -o OUTDIR, --outdir OUTDIR
                        name of output directory (will be created if does not exist) (default: output_pesci/)

  --label_species1 LABEL_SPECIES1
                        name of species1 (label for plots & outputs) (default: sp1)

  --label_species2 LABEL_SPECIES2
                        name of species2 (label for plots & outputs) (default: sp2)

  -f FIGURE_FORMAT, --figure_format FIGURE_FORMAT
                        {svg,png,pdf} format for output figures (default: svg)

  --plot_thresh PLOT_THRESH
                        threshold to plot matches in resulting comparison matrix & search for co-expressed marker gene pairs (default: 0)

  --min_fc MIN_FC       minimum fold-change to be considered marker of a cluster - only used for the co-
                        expressed marker gene table (to set in conjunction with --marker_specificity)
                        (default: 1.5)