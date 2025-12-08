Required input files
=====================

Important in data preperation to have same gene id in the orthology and matrix files AND that gene ids are unique for each species (recommend adding prefix with sp name). go to how to run for command line arguments to use.

Orthology file
--------------
Same for all

Single-cell expression data
----------------------------

Different accepted input formats, and can be different for each of the two study species

1 - Sparse count matrix and cell to cluster annotation
-------------------------------------------------------

For instance, one way to format from a Seurat object:

```
library(Matrix)
library(R.utils)
library(data.table)
library(tidyverse)

data_dir <- 'scdata_folder/'

Preferred: 
counts <- GetAssayData(mySeuratObj, assay="RNA", slot='counts')

OR

counts <- mySeuratObj@assays$RNA@counts

OR 

counts <- mySeuratObj[["RNA"]]$counts

# Write counts matrix
writeMM(counts, paste0(data_dir, 'matrix.mtx'))
gzip(paste0(data_dir, 'matrix.mtx'))

# Write cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes), paste0(data_dir, 'barcodes.tsv'),
           col_names = FALSE)
gzip(paste0(data_dir, 'barcodes.tsv'))

# Write feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", paste0(data_dir, 'features.tsv'),
           col_names = FALSE)
gzip(paste0(data_dir, 'features.tsv'))

# Write cell to cluster
write.table(mySeuratObj$Cluster_labels, file = "/path/to/Cluster.tsv", sep = "\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

```

2 - Scanpy H5ad file
----------------------

same as 1. but after running 1 do an extra-step in python:

import scanpy as sc
import pandas as pd
import numpy as np

fwdata = sc.read('Lvar.counts.mtx').transpose()
cellnames = pd.read_csv('Lvar.col.txt', header=None, index_col=0)
genenames = pd.read_csv('Lvar.row.txt', header=None, index_col=0)
cellnames.index.name = 'index'
genenames.index.name = 'index'
fwdata.obs = cellnames
fwdata.var = genenames

cluster_assignments = pd.read_csv('Lvar.tsv', header=None, sep='\t', index_col=0, names=['cluster id'])
cluster_assignments.index.name = 'index'

fwdata.obs=fwdata.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)

fwdata.write_h5ad('Lvar.h5ad')


3 - Dense count matrix and cell to cluster annotation
-------------------------------------------------------








Optional input files
=====================