# Help With Input Files

## Required Input Files

Pesci requires 3 main inputs: the **single-cell gene expression data** for each of the two species under comparison and a **gene orthology file** listing all orthologous gene pairs (1-to-1 and many-to-many) across these two species. Single-cell gene expression data can be provided in any of the following formats: (1) a **sparse expression matrix** (CellRanger-like directory), (2) a **scanpy h5ad** file **OR** (3) a **dense expression matrix**. On this page, we provide example code showing how to convert a **Seurat object** into one of these accepted formats. Depending on input format, an additional file giving the **cell barcode to cell cluster correspondence** might be necessary.

Also say that if several lib and one file per lib --> possible to give them all, see examples page.

> [!WARNING]
> Two important aspects of data preparation for pesci are to ensure (i) that **gene ids** (or gene names) are the **same across the provided gene expression matrix and gene orthology file** and (ii) that gene ids are **unique to each species** (if unsure, we recommend adding a species prefix to gene names, i.e. for instance Procro_TTN and Cragig_TTN).

### Orthology file

Same for all
Tab or comma delimiter gz or not 
Could be from biomart 
Could be from broccoli or orthofinder (give file names)

### Single-cell expression data

In this section, we describe the different accepted single-cell expression data formats and how they can be generated from a Seurat Object.

> [!TIP]
> To the exception of scanpy's .h5ad, all input files (single-cell data, orthology file, cell-to-cluster files) can be compressed in .gz or .gz2 or left uncompressed. Similarly, where applicable (i.e. not h5ad or sparse single-cell data), both tab and comma-separated formats are accepted.

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


Add code to add prefix species (if in orthology file) 

gene_names <- paste("prefix_", gene_names, sep = "")

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

Alternatively, for datasets processed with scanpy, pesci can directly work with h5ad (but needs the raw counts, if no layer 'counts' will check the data.X but will error if not integers)

```python
import scanpy as sc
import pandas as pd

data = sc.read_10x_mtx('Cragig_sparse_matrix/')
data.layers['counts'] = data.X.copy()

cluster_assignments = pd.read_csv('data/Cragig_cell_id.tsv', sep='\t', index_col=0, usecols=[0, 1])
data.obs = data.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)

data.write_h5ad('data/Cragig_matrix.h5ad')
```

Add code to add prefix species (if in orthology file) and/or to substitute gene names using a table


**3 - Dense count matrix and cell to cluster annotation table**
no code, less optimal (if Seurat object recommend making a sparse matrix or h5ad). Only supported to easily load GEO datasets.
import to use.csv is comma-sep and .tsv for tab-sep

Add code to add prefix species (if in orthology file) and/or to substitute gene names using a table


## Optional Input Files

Optionally, 
Broad annotations (additional column in the cluster file or h5ad)
