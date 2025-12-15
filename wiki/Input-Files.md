# Help With Input Files

## Required Input Files

Pesci requires 3 main inputs: a **gene orthology file** listing all orthologous gene pairs (1-to-1 and many-to-many) between the two species under comparison and the **single-cell gene expression count data** for each of these two species.

Single-cell gene expression counts can be provided in any of the following formats: (1) a **sparse expression matrix** (a CellRanger-like directory), (2) a **scanpy h5ad** file **OR** (3) a **dense expression matrix**. It is also possible to provide several input files for a single species, for instance in cases where several libraries were sequenced and stored in distinct files (please see example run in [Pesci Examples](https://github.com/eparey/pesci/blob/main/wiki/Examples.md)). Depending on input format, an additional file giving the **cell barcode to cell cluster correspondence** might be necessary.

On this page, we describe all accepted input formats and provide example code to convert **Seurat Objects** into these formats. 


> [!WARNING]
> Two important aspects of data preparation for pesci are to ensure (i) that **gene ids** (or gene names) are the **same across the provided gene expression matrix and gene orthology file** and (ii) that these gene ids are **unique to each species** (if unsure, we recommend adding a species prefix to gene names, *i.e.* for instance Procro_TTN and Cragig_TTN for the gene encoding TTN).

### Orthology file

The gene orthology file is a **two-columns** file, either **tab- (.tsv)** or **comma-separated (.csv)**. It can be compressed (.gz or .gz2) or uncompressed. Apart from these two rules, the expected format is flexible, accommodating inputs from different sources.

The easiest way to generate this file is to use pre-computed orthologies, provided the two species under comparisons are available in comparative genomics databases (for instance [Ensembl](www.ensembl.org/)). Alternatively, tools like [OrthoFinder](https://github.com/davidemms/OrthoFinder) and [Broccoli](https://github.com/rderelle/Broccoli) can produce orthologs for user-specific datasets.

Pesci accepts "simple" .csv or .tsv files, .tsv generated from [Ensembl bioMart](www.ensembl.org/info/data/biomart/index.html), files in [Broccoli](https://github.com/rderelle/Broccoli) format, or in [OrthoFinder](https://github.com/davidemms/OrthoFinder) format (more details in examples below). A real example (Oyster-larvae vs Flatworm-larvae comparison, in [Broccoli](https://github.com/rderelle/Broccoli)) is provided in the [data folder](https://github.com/eparey/pesci/blob/main/data/orthologous_pairs_Procro-Cragig.txt).

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

- **Example 2:**  [Ensembl bioMart](www.ensembl.org/info/data/biomart/index.html) (or [Ensembl Metazoa bioMart](www.metazoa.ensembl.org/info/data/biomart/index.html)
    
    (two species are in ensembl or ensembl metazoa or any other ensembl collection)
    How to get = biomart, select species 1 Dataset click attributes tick the circle Homologues (Max select 6 orthologues), in the gene menu select only 1 type of id Gene stable ID (or protein if matches single cell) + Orthologue select species 2 only gene id (or any) + click results and tsv and unique results only)

    Tab separated + possible empty col2 (can be removed using biomart filter but pesci does not care)

    Ensembl bioMart format (tab-separated)
    ```
    Procro_TTN  Cragig_TTN
    Procro_g1332    Cragig_g145
    Procro_g1332    Cragig_g146
    Procro_g1333    
    Procro_g1334    
    Procro_g1335    Cragig_g5587
    Procro_g1336    Cragig_g5587
    ```

- **Example 3:** [Broccoli](https://github.com/rderelle/Broccoli) orthologous gene pairs file (~ `BroccoliOUTdir_step4/orthologous_gene_pairs.txt`):
    Broccoli is a tool to infer orthologous gene groups and orthologous gene pairs using a mixed phylogeny-network approach. To accommodate Broccoli inputs, pesci allows for genes from other species to be present in the orthology file (i.e. not just the two species under comparisons). In addition, genes from the same species do not have to stick to being in the same column throughout the file.

    Broccoli format (tab-separated)
    ```
    Cragig_TTN  Procro_TTN
    Procro_g1332    Cragig_g145
    Procro_g1332    Cragig_g146
    Pecmax_TTN1 Procro_TTN
    Procro_g1335    Cragig_g5587
    Cragig_g5587    Procro_g1336    
    ```


- **Example 4:** [OrthoFinder](https://github.com/davidemms/OrthoFinder) (~ `OrthoFinderOUT/Orthologues/Species1__v__Species2.csv`):


    [OrthoFinder](https://github.com/davidemms/OrthoFinder) finds orthogroups and orthologous gene pairs using phylogenetic gene tree reconstruction. 
    need to remove the first column in the orthofinder orthologue file:
    `cut -f 2,3 Species1__v__Species2.csv > orthofinder_orthologues_sp1-sp2_ok.tsv`

    Tab separated + single gene or a comma-separated list of genes (from the same species!) can be present in col2
    ```
    Procro_TTN  Cragig_TTN
    Procro_g1332    Cragig_g145, Cragig_g146
    Procro_g1335, Procro_g1336    Cragig_g5587
    ```

> [!IMPORTANT]
> Any combination of these formats will be accepted (for instance, in-file species columns swap and/or additional species present in an OrthoFinder-like file). The only requirements are: the file must have exactly two columns, gene ids must be the same in the orthology file and in the single cell count matrices and gene ids must be unique to each species (and unique with respect to genes of other species potentially also present in the file).

### Single-cell expression data

In this section, we describe the different accepted single-cell expression data formats and how they can be generated from Seurat Objects.

> [!TIP]
> To the exception of scanpy's .h5ad, all input files (single-cell data, orthology file, cell-to-cluster files) can be compressed in .gz or .gz2 or left uncompressed. Similarly, where applicable (i.e. not h5ad or sparse single-cell data), both tab (.tsv) and comma-separated (.csv) formats are accepted.

**1. Sparse count matrix and cell to cluster annotation table**

Describe the format
Directory three files: matrix.mtx barcodes.tsv features.tsv
Cell to cluster separator can be tab or comma (automatically detected) + second column or any name provided with colclust
--> can also retain only a subset of clusters filter_out and keep_only args (see options (Inputs arg) and examples)

Explain (would it remove cells not in clusters?, i.e. if giving the CellRanger dir + clusters --> yes but I need to check performance)

The following R code is provided as an example to format a Seurat Object into a Sparse Matrix for pesci. It creates a CellRanger-like directory (hereafter named "Cragig_sparse_matrix/") that can be directly provided as argument to --matrix1 or --matrix2.

```R
library(Matrix)
library(R.utils)
library(data.table)
library(Seurat)

# Load the Seurat object
mySeuratObj <- readRDS("Cragig_seurat_object.rds")

# Set and create the output directory (to store the sparse matrix files)
data_dir <- 'Cragig_sparse_matrix/'
dir.create(data_dir)

# Get the counts matrix from the Seurat object
counts <- GetAssayData(mySeuratObj, assay="RNA", slot='counts') #Seurat v4 and v5
#counts <- mySeuratObj@assays$RNA@counts #Seurat v3 and v4 only
#counts <- mySeuratObj[["RNA"]]$counts  #Seurat v5 only

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

# Write the cell-to-cluster annotation (replacing 'cluster_labels' by the name of the meta.data column containing clusters info if different, e.g 'seurat_clusters')
write.table(mySeuratObj$cluster_labels, file="Cragig_cell_clusters.tsv", sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
```

#optional: add a species prefix to gene names



Add code to add prefix species (if in orthology file) 

gene_names <- paste("prefix_", gene_names, sep = "")


#optional: substitute gene names using a conversion table
and/or to substitute gene names using a table

> library(hashmap)
read.table stringsAsFactors = FALSE
> keys = c("a", "b", "c", "d")
> values = c("A", "B", "C", "D")
> my_dict <- hashmap(keys, values)
> my_vect <- c("b", "c", "c") 
add check that all elements are in the dict
> translated <- my_dict$find(my_vect)
> translated
[1] "B" "C" "C"


**2. Scanpy h5ad file**

Alternatively, for datasets processed with scanpy, pesci can directly work with h5ad files, provided the raw counts are stored in the object.  (but needs the raw counts, if no layer 'counts' will check the data.X but will throw InputError if not integers)

The following code snippet shows how to generate and h5ad from a CellRanger(-like) directory (as generated in the previous section (link) or as CellRanger). Not-necessary if already done 1, just an alternative!!

```python
import scanpy as sc
import pandas as pd

#load CellRanger(-like) directory in scanpy
data = sc.read_10x_mtx('Cragig_sparse_matrix/')
data.layers['counts'] = data.X.copy() #in this case we know for sure that data.X are raw counts, so we store it in the count layer

#load cluster annotations and add to the object
cluster_assignments = pd.read_csv('data/Cragig_cell_id.tsv', sep='\t', index_col=0, usecols=[0, 1])
data.obs = data.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)

#optional: add a species prefix to gene names

#optional: substitute gene names using a conversion table

data.write_h5ad('data/Cragig_matrix.h5ad')
```



**3 - Dense count matrix and cell to cluster annotation table**
no code, less optimal (if Seurat object recommend making a sparse matrix or h5ad). Only supported to easily load GEO datasets.
import to use.csv is comma-sep and .tsv for tab-sep

Add code to add prefix species (if in orthology file) and/or to substitute gene names using a table
 same as 1 Cell to cluster separator can be tab or comma (automatically detected) + second column or any name provided with colclust
--> can also retain only a subset of clusters filter_out and keep_only args (see options (Inputs arg) and examples)


## Optional Input Files

Optionally, 
Broad annotations (additional column in the cluster file or h5ad)


## References

Broccoli

OrthoFinder

Ensembl BioMart