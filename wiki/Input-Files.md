# Help With Input Files

## Required Input Files

Important in data preperation to have same gene id in the orthology and matrix files AND that gene ids are unique for each species (recommend adding prefix with sp name). go to how to run for command line arguments to use.

### Orthology file

Same for all

> [!WARNING]
> Potentially most tricky = get gene names the same in matrix and orthology (maybe give a few pointer codes)

### Single-cell expression data

Different accepted input formats, and can be different for each of the two study species

> [!TIP]
> All input files, to the exception of scanpy's .h5ad files, can be compressed in .gz and will be properly read by pesci (whether compressed or not).

**1. Sparse count matrix and cell to cluster annotation table**

Explain (would it remove cells not in clusters?, i.e. if giving the cellrnager dir)
The following R code is provided as an example to format a Seurat Object into a Sparse Matrix for pesci. It creates a CellRanger-like directory (hereafter named "Cragig_sparse_matrix/") that can be directly provided as argument to --matrix1 or --matrix2 in pesci.

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
counts <- GetAssayData(mySeuratObj, assay="RNA", slot='counts') #seurat v4 and v5
#counts <- mySeuratObj@assays$RNA@counts #seurat v3 and v4 only
#counts <- mySeuratObj[["RNA"]]$counts  #seurat v5 only

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
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", paste0(data_dir, 'features.tsv'),
           col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))

# Write the cell-to-cluster annotation (replacing 'cluster_labels' by the name of the meta.data column containing clusters info if different, e.g 'seurat_clusters')
write.table(mySeuratObj$cluster_labels, file = "Cragig_cell_clusters.tsv", sep = "\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
```

**2. Scanpy h5ad file**

Alternatively, for datasets processed with scanpy, pesci can 

```python
import scanpy as sc
import pandas as pd

data = sc.read_10x_mtx('Cragig_sparse_matrix/')
data.layers['counts'] = data.X.copy()

cluster_assignments = pd.read_csv('data/Cragig_cell_id.tsv', sep='\t', index_col=0, usecols=[0, 1])
data.obs = data.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)

data.write_h5ad('data/Cragig_matrix.h5ad')
```


**3 - Dense count matrix and cell to cluster annotation table**
no code, less optimal (if seurat object recommend making a sparse matrix). Only supported to easily load GEO datasets.



## Optional Input Files

Broad annotations (additional column in the cluster file or h5ad)
