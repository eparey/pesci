pesci - Pretty Easy Single-cell Comparisons using ICC
======================================================

.. image:: https://img.shields.io/badge/python-3.7+-blue.svg
    :target: https://www.python.org/downloads/

Pesci is an efficient and user-friendly implementation of the Iterative Comparison of Coexpression (ICC) algorithm applied to the **comparison of single-cell datasets across two study species** (`Najle, Grau-Bové et al. <https://doi.org/10.1016%2Fj.cell.2023.08.027>`_).

Pesci takes as input **single-cell expression count matrices** (raw count matrices, CellRanger directory and/or h5ad files), **cell cluster annotations** and **gene orthologies** files.

To calculate expression similarity, pesci uses 1-to-1 orthologs and accounts for  many-to-many orthologs by selecting the gene pair with the most similar expression (i.e. most likely to reflect the ancestral gene function).

.. toctree::
   :caption: Manual
   :name: manual
   :maxdepth: 1

   getting_started.rst
   input_description.rst
   output_advanced.rst


.. toctree::
   :caption: Project Information
   :name: project_information
   :maxdepth: 1

   modules.rst
