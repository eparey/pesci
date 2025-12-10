List of available options
==========================

This page provides a description of all available options when running pesci. For more details on input and output files, including R code to prepare input files from Seurat objects, please see `Help With Input Files <https://github.com/eparey/pesci/blob/main/wiki/Input-Files.md>`_ and []. For example use cases with specific options, please see [].


Minimal command-line examples
-----------------------------

Below are minimal examples to run pesci, each showcasing different accepted formats for input single-cell expression datasets. The corresponding data files are available as examples in the data/ directory.

  - **Example 1**: species 1 data as a sparse count matrix (same format as a cellranger directory) and species 2 data as a scanpy .h5ad

  .. code-block:: sh

    pesci --matrix1 data/Cragig_sparse/ --matrix2 data/Procro_matrix_EM.h5ad --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt

  - **Example 2**: data for both species provided as dense count matrices (least optimal format)

  .. code-block:: sh

    pesci --matrix1 data/Cragig_matrix_EM.tsv.gz --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt

  - **Example 3**: data for both species provided as dense count matrices, but with species 2 data from 3 different libraries with counts respectively saved in 3 different files

  .. code-block:: sh

    pesci --matrix1 data/Cg_matrix_EM_part1.tsv data/Cg_matrix_EM_part2.tsv data/Cg_matrix_EM_part3.tsv --matrix2 data/Procro_matrix_EM.tsv.gz --clusters1 data/Cragig_cell_id.tsv --clusters2 data/Procro_cell_id.tsv -g data/orthologous_pairs_Procro-Cragig.txt

*Note*: files can be .gz or not. Any combination of valid input formats for species 1 and 2 is accepted.

*important* note somewhere: if running different comparisons use sp. labels and/or different output directory, otherwise with warning will reuse saying already computed


Usage
------

All available options are described below and can be printed using `pesci --help`. Note that flags (--force, --ono2one_only, --do_not_plot_warn, --show_auto_threshold, --no_pbar) do not require any arguments to be set, for instance `` will effectively force recomputation of all intermediary files.

**options:**
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

**required arguments:**
  --matrix1 MATRIX1
                        gene expression per cell for species 1 (raw count matrix), either: a dense matrix text file (can be .gz), a cellranger directory or an h5ad file; several inputs can be specified if multiple libraries

  --matrix2 MATRIX2
  						          gene expression per cell for species 2 (raw count matrix), either: a dense matrix text file (can be .gz), a cellranger directory or an h5ad file; several inputs can be specified if multiple libraries

  --clusters1 CLUSTERS1
                        cell-to-clusters table for species 1 or name of cluster column in h5ad

  --clusters2 CLUSTERS2
                        cell-to-clusters table for species 2 or name of cluster column in h5ad

  -g ORTHOLOGOUS_GENE_PAIRS, --orthologous_gene_pairs ORTHOLOGOUS_GENE_PAIRS
                        gene orthologies file (tab-delimited, 1 pair per line)

**recommended arguments:**
  -o OUTDIR, --outdir OUTDIR
                        name of output directory (will be created if does not exist) (default: output_pesci/)

  --label_species1 LABEL_SPECIES1
                        name of species1 (label for plots & outputs) (default: sp1)

  --label_species2 LABEL_SPECIES2
                        name of species2 (label for plots & outputs) (default: sp2)

**resources:**
  -c CORES, --cores CORES
                        number of cores to use for parallel operations (default: 1)

**running:**
  --force               
                        recompute all intermediary files (per cluster normalized gene expression and ec scores) even if outputs already exist (default: False)

  --seed SEED           
                        set random seed (default: 123)

  --marker_specificity MARKER_SPECIFICITY
                        marker specificity (btwn 0.5 and 0.95) - used in fold change expression calculation
                        as follows: 0.5 computes fold change in cluster 'a' vs median expression across
                        clusters, 0.75 computes fold change in cluster 'a' vs 75 expression percentile across
                        clusters (3rd quartile) (default: 0.5)

  --ono2one_only        
                        use 1-to-1 orthologs only (default: False)

  --ec_threshold_many EC_THRESHOLD_MANY
                        if set, instead of using only the best pair for many-to-many and many-to-one orthologs, retain all distinct pairs with score > threshold
                        (default: None)

  --do_not_downsample
                        use all 1-to-1 orthologs (instead of 1000 randomly chosen) for selection of best pairs for the
                        many-to-many and many-to-one orthologs (slower) (default: False)

  --random_id RANDOM_ID
                        triggers a randomization of orthologies to run pesci on randomized orthologs; stores results in outdir/random_id
                        (default: )
  --no_pbar             
                        do not show progress bar (default: False)


**inputs:**
  --min_umi MIN_UMI     
                        minimum total umi for a gene to be retained for comparison (default: 10)
  
  --colclust COLCLUST   
                        name of the column corresponding to clusters in the cell to cluster annotation file
                        (does not apply to h5ad input as it is specified by --clusters1/--cluster2 directly)
                        - if not set, will use column `cluster_name`, if it does not exist, will attempt to
                        use the 2nd column

                        use --colclust1 or --colclust2 to specify different column names
                        for species1 or species2. (default: cluster_name)
  
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

                        use --colbroad1 or --colbroad2 to specify different column names
                        for species1 or species2. (default: None)
  
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

**outputs:**
  -f FIGURE_FORMAT, --figure_format FIGURE_FORMAT
                        {svg,png,pdf} format for output figures (default: pdf)

  -r REORDER, --reorder REORDER
                        method for ordering cell clusters on the heatmap: DiagKeep (default): try to maximise matches on the diagonal, keeping input order as much as possible; Clust: use hierarchical clustering of rows and columns (average linkage, euclidean distance) or None: keep input order. This option is ignored if providing broad annotations. (default: DiagKeep)

  --min_fc MIN_FC       
                        minimum fold-change to be considered marker of a cluster - only used for the co-
                        expressed marker genes table (to set in conjunction with --marker_specificity).
                        (default: 1.5)

  --seaborn_cmap SEABORN_CMAP
                        name of the seaborn colormap for the heatmap. (default: BuPu)

  --show_auto_threshold 
                        experimental option to highlight matches above thresholds on the heatmap plot. (default: False)

  --do_not_plot_warn
                        do not plot warning on the heatmap for matches driven by < 10 genes. (default: True)