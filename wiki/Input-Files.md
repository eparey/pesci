# Help With Input Files

This page describes accepted input file formats and how to produce them.

## Table of Content

- [Required Input Files](#required-input-files)
    - [Orthology file](#orthology-file)
    - [Single-cell expression counts](#single-cell-expression-counts)
- [Optional Input Files](#optional-input-files)
- [References](#references)


## Required Input Files

Pesci requires 3 main inputs: a **gene orthology file** listing all orthologous gene pairs (one-to-one and many-to-many) between the two species under comparison and the **single-cell gene expression count data** for each of these two species.

Single-cell gene expression counts can be provided in any of the following formats: a **sparse expression matrix** (a CellRanger-like directory), a **scanpy h5ad** file **OR**  a **dense expression matrix**. Several input files for a single species can be provided, for instance in cases where several libraries were sequenced and stored in distinct files (please see an example run in [Pesci Examples](https://github.com/eparey/pesci/blob/main/wiki/Examples.md)). Depending on input format, an additional file giving the **cell barcode to cell cluster correspondence** might be necessary.

On this page, we describe all accepted input formats and present example code to convert **Seurat Objects** into these formats. 


> [!IMPORTANT]
> Two important aspects of data preparation for pesci are to ensure (i) that **gene ids** (or gene names) are the **same across the provided gene expression matrix and gene orthology file** and (ii) that these gene ids are **unique to each species** (if unsure, we recommend adding a species prefix to gene names, *i.e.* for instance Procro_TTN and Cragig_TTN for the gene encoding TTN, see below for formatting help).

### Orthology file

The gene orthology file is a **two-columns** file, either **tab- (.tsv)** or **comma-separated (.csv)**. It can be compressed (.gz or .gz2) or uncompressed. Apart from these two rules, the expected format is flexible, accommodating inputs from different sources.

The easiest way to generate this file is to use pre-computed orthologies, provided the two species under comparisons are available in existing comparative genomics databases (for instance [Ensembl](www.ensembl.org/)). Alternatively, tools like [OrthoFinder](https://github.com/davidemms/OrthoFinder) and [Broccoli](https://github.com/rderelle/Broccoli) can infer orthologs on user-specific datasets.

Pesci accepts: "simple" .csv or .tsv files, .tsv generated from [Ensembl bioMart](www.ensembl.org/info/data/biomart/index.html), files in [Broccoli](https://github.com/rderelle/Broccoli) format or in [OrthoFinder](https://github.com/davidemms/OrthoFinder) format. Please find more details in examples below:


- **Example 1:** simple .csv or .tsv

    csv
    ```
    Procro_TTN,Cragig_TTN
    Procro_g1332,Cragig_g145
    Procro_g1332,Cragig_g146
    Procro_g1335,Cragig_g5587
    Procro_g1336,Cragig_g5587
    ```

    tsv
    ```
    Procro_TTN  Cragig_TTN
    Procro_g1332    Cragig_g145
    Procro_g1332    Cragig_g146
    Procro_g1335    Cragig_g5587
    Procro_g1336    Cragig_g5587
    ```

- **Example 2:**  [Ensembl bioMart](www.ensembl.org/info/data/biomart/index.html) (or [Ensembl Metazoa bioMart](www.metazoa.ensembl.org/info/data/biomart/index.html))
    
    Orthologous genes for a pair of species can be obtained from Ensembl bioMart as follows: 1. choose dataset: select genes of species 1, 2. click attributes and tick the "Homologues (Max select 6 orthologues)" circle, 3. expand the gene menu and select only a single type of ID (corresponding to IDs in the single-cell data), 4. click the orthologues for species 2, selecting again only one type of IDsm 5. click "Results" on the top left of the screen, tick unique results only and click Go to get the orthologies in tsv format.

    Ensembl bioMart format (tab-separated, possibly empty column 2)
    ```
    Procro_TTN  Cragig_TTN
    Procro_g1332    Cragig_g145
    Procro_g1332    Cragig_g146
    Procro_g1333    
    Procro_g1334    
    Procro_g1335    Cragig_g5587
    Procro_g1336    Cragig_g5587
    ```

- **Example 3:** [Broccoli](https://github.com/rderelle/Broccoli) orthologous gene pairs file (`BroccoliOUT/dir_step4/orthologous_gene_pairs.txt`):
    
    Broccoli is a tool to infer orthologous gene groups and orthologous gene pairs using a mixed phylogeny-network approach. To accommodate Broccoli inputs, pesci allows for genes from other species to be present in the orthology file (i.e. not just the two species under comparisons). In addition, genes from the same species do not have to stick to being in the same column throughout the file.

    Broccoli format (tab-separated, additional species + column swap)
    ```
    Cragig_TTN  Procro_TTN
    Procro_g1332    Cragig_g145
    Procro_g1332    Cragig_g146
    Pecmax_TTN1 Procro_TTN
    Procro_g1335    Cragig_g5587
    Cragig_g5587    Procro_g1336    
    ```

    A real example (Oyster-larvae vs Flatworm-larvae comparison, in [Broccoli](https://github.com/rderelle/Broccoli) format) is provided in the [data folder](https://github.com/eparey/pesci/blob/main/data/orthologous_pairs_Procro-Cragig.txt).


- **Example 4:** [OrthoFinder](https://github.com/davidemms/OrthoFinder) (`OrthoFinderOUT/Orthologues/Species1__v__Species2.csv`):


    [OrthoFinder](https://github.com/davidemms/OrthoFinder) builds orthology groups and finds orthologous gene pairs using phylogenetic gene tree reconstruction. 
    
    To be used with pesci the first column must be removed from the OrthoFinder orthologous pairs file:
    `cut -f 2,3 Species1__v__Species2.csv > orthofinder_orthologues_sp1-sp2_ok.tsv`

    OrthoFinder format (tab-separated, single gene or a comma-separated list of genes in each column)
    ```
    Procro_TTN  Cragig_TTN
    Procro_g1332    Cragig_g145, Cragig_g146
    Procro_g1335, Procro_g1336    Cragig_g5587
    ```

> [!IMPORTANT]
> Any combination of these formats will be accepted (for instance, in-file species columns swap and/or additional species present in an OrthoFinder-like file). The only requirements are: the file must have exactly two columns, gene ids must be the same in the orthology file and in the single cell count matrices and gene ids must be unique to each species (and unique with respect to genes of other species potentially also present in the file).

### Single-cell expression counts

In this section, we describe the different accepted single-cell expression data formats and how they can be generated from Seurat Objects. Note that several input files for a single species can be provided, for cases where several libraries were sequenced and stored in distinct files (please see an example run in [Pesci Examples](https://github.com/eparey/pesci/blob/main/wiki/Examples.md)).

> [!TIP]
> To the exception of scanpy's .h5ad, all input files (single-cell data, orthology file, cell-to-cluster files) can be compressed in .gz or .gz2 or left uncompressed. Similarly, where applicable (i.e. not .h5ad or sparse single-cell data), both tab- (.tsv) and comma-separated (.csv) formats are accepted.

- **1. Sparse count matrix and cell to cluster annotation table**

    Pesci accepts sparse count matrices formatted as a directory containing: the count matrix (matrix.mtx), the barcodes i.e. column names (barcodes.tsv) and the gene names i.e. row names (features.tsv). This is the exact same format as a CellRanger output directory, and can also be easily generated from a Seurat Object (see below). The directory can be provided as argument to the --matrix1 (or --matrix2 for species 2) argument.

    In addition, a file giving the cell barcode to cluster correspondence is also required (--clusters1 or --clusters2 argument). This file can either be a comma-separated (.csv) or tab-separated (.tsv) file, with any number of columns, with cell barcodes in the first column and cluster annotation in any of the other columns. By default, pesci uses the second column as cluster annotation, but this behaviour can be overruled by providing a column name to the option --colclust1 (or --colclust2). Note that any cell barcode that is not found in the cell-to-cluster file will be ignored. Additional arguments can be specified in order to use only a subset of the cluster annotations, if necessary (see --filter_out and --keep_only in [Pesci Usage](https://github.com/eparey/pesci/blob/main/wiki/How-to-run.md) and [Pesci Examples](https://github.com/eparey/pesci/blob/main/wiki/Examples.md)). 

    The following R code shows how to format a Seurat Object into a sparse matrix for pesci. It creates a CellRanger-like directory (hereafter named "Cragig_sparse_matrix/") that can be directly provided as argument to --matrix1 (or --matrix2) along with the corresponding barcode-to-cluster file for --clusters1 (or --clusters2).

    ```R
    library(Matrix)
    library(R.utils)
    library(data.table)
    library(Seurat)
    library(readr)

    # Load the Seurat object
    mySeuratObj <- readRDS("Cragig_seurat_object.rds")

    # Set and create the output directory (to store the sparse matrix files)
    data_dir <- 'Cragig_sparse_matrix/'
    dir.create(data_dir)

    # Get the counts matrix from the Seurat object
    counts <- GetAssayData(mySeuratObj, assay="RNA", slot='counts') #Seurat v4 and v5
    #counts <- mySeuratObj@assays$RNA@counts #Seurat v3 and v4 only
    #counts <- mySeuratObj[["RNA"]]$counts  #Seurat v5 only, may need to JoinLayers before if several 

    # Write the counts matrix
    writeMM(counts, paste0(data_dir, 'matrix.mtx'))
    gzip(paste0(data_dir, 'matrix.mtx'))

    # Write the cell barcodes
    barcodes <- colnames(counts)
    write_delim(as.data.frame(barcodes), paste0(data_dir, 'barcodes.tsv'),
               col_names = FALSE)
    gzip(paste0(data_dir, 'barcodes.tsv'))

    # Write the feature names
    gene_names <- rownames(counts)
    features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
    write_delim(as.data.frame(features),delim = "\t", paste0(data_dir, 'features.tsv'),
               col_names = FALSE)
    gzip(paste0(data_dir, 'features.tsv'))

    # Write the cell-to-cluster annotation 
    #(replacing 'cluster_labels' by the name of the meta.data column containing clusters info if different, e.g 'seurat_clusters')
    write.table(mySeuratObj$cluster_labels, file="Cragig_cell_clusters.tsv", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
    ```

    - Optional alternative: renaming genes with a species prefix if necessary for uniqueness and/or consistency with the orthology file:

        ```R
        #optional: add a species prefix to gene names
        gene_names <- rownames(counts)
        gene_names <- paste("Cragig_", gene_names, sep = "")
        features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
        write_delim(as.data.frame(features),delim = "\t", paste0(data_dir, 'features.tsv'),
                   col_names = FALSE)
        gzip(paste0(data_dir, 'features.tsv'))
        ```

    - Optional alternative 2: renaming genes using a conversion table if necessary for uniqueness and/or consistency with the orthology file:

        ```R
        library(hashmap)

        gene_names <- rownames(counts)

        #optional: substitute gene names using a conversion table
        read.table stringsAsFactors = FALSE
        keys = c("a", "b", "c", "d")
        values = c("A", "B", "C", "D")
        my_dict <- hashmap(keys, values)
        my_vect <- c("b", "c", "c")

        translated <- my_dict$find(my_vect)
        translated
        [1] "B" "C" "C"
        ```

- **2. Scanpy h5ad file**

    Alternatively, for datasets processed with scanpy, pesci can directly work with .h5ad files, provided the raw counts are stored in the object. In practice, pesci will first check for the presence of a 'counts' layer, then for data.raw.X, and if none of these are found in the object pesci will check that the data.X contains integer count values, if yes pesci will use it, otherwise an error will be returned.

    Cluster annotations are expected to be stored in the scanpy AnnData object saved as .h5ad, so no additional file is necessary - simply specify the name of the annotation column within the object to --clusters1 (or --clusters2). 

    The following code snippet shows how to generate an .h5ad file from a CellRanger(-like) directory (as generated in the previous section or as produced by CellRanger), as an alternative accepted format.

    ```python
    import scanpy as sc
    import pandas as pd

    #load CellRanger(-like) directory in scanpy
    data = sc.read_10x_mtx('Cragig_sparse_matrix/')
    data.layers['counts'] = data.X.copy() #in this case we know that data.X are the raw counts, we store it in the counts layer

    #load cluster annotations and add to the AnnData object
    cluster_assignments = pd.read_csv('data/Cragig_cell_id.tsv', sep='\t', index_col=0, usecols=[0, 1])
    data.obs = data.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)
    #print(data.obs) #check newly-added column 'cluster_name'

    #optional: add a species prefix to gene names
    data.var.index = 'Cragig_' + data.var.index.astype(str)

    #optional: substitute gene names using a conversion table

    data.write_h5ad('data/Cragig_matrix.h5ad')
    ```



- **3. Dense count matrix and cell to cluster annotation table**

    This format is the least optimal format for pesci, but allows to more easily use raw dense matrices as found in GEO datasets. The matrix file can be tab- (.tsv) or comma (.csv) separated, with genes in rows and cell barcodes in columns (first column contains the gene names and first row the cell barcodes). A real example is shown in the [data folder](https://github.com/eparey/pesci/blob/main/data/Cragig_matrix_EM.tsv.gz). An additional cell to cluster annotation file is also required, see 1. on sparse matrix above for details.

    Add code to add prefix species (if in orthology file) and/or to substitute gene names using a table

    Such a count matrix **can be provided to pesci as is**. The code below is only useful for cases where it is necessary to rename genes with a species prefix or a conversion table. The code requires pesci to be installed and saves the matrix with corrected gene names as an .h5ad file (more optimal with pesci).

    ```python
    import scanpy as sc
    import pandas as pd
    from pesci import normalize as pnm

    expr_mat = 'data/Cragig_matrix_EM.tsv' #dense matrix file
    cores = 4 #number of threads to use for loading
    fmt = 'tsv' #change to csv if file is comma-separated

    expr, genes, cells = pnm.load_matrix_dense(expr_mat, cores, fmt) #if the file is .gz, add the argument open_func=gzip.open

    #Optional: add a species prefix to gene names
    genes = ['Cragig_'+i for i in genes]

    #OR use a conversion table to rename genes
    #

    data = sc.AnnData(expr.T, pd.DataFrame(index=cells), pd.DataFrame(index=genes))

    clusters = pd.read_csv('data/Cragig_cell_id.tsv', sep='\t', index_col=0, usecols=[0, 1])
    data.obs = data.obs.merge(clusters, how='left', left_index=True, right_index=True)

    data.write_h5ad('data/Cragig_matrix_raw.h5ad')
    ```

## Optional Input Files

Optionally, broad annotations can be provided on top of cluster annotations and are used to order the output single-cell comparison heatmap. These can be given as an additional column in the cell-to-cluster file or in the AnnData object saved as .h5ad (in both case, the name of the column is to be specified with --colbroad).


## References

- Broccoli

- OrthoFinder

- Ensembl BioMart

- Seurat 

- Scanpy