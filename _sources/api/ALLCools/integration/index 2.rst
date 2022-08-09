:py:mod:`ALLCools.integration`
==============================

.. py:module:: ALLCools.integration


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cca/index.rst
   confusion/index.rst
   harmony/index.rst
   seurat/index.rst


Package Contents
----------------

.. py:function:: simple_cca(adata, group_col, n_components=50, random_state=0)


.. py:function:: incremental_cca(a, b, max_chunk_size=10000, random_state=0)

   Perform Incremental CCA by chunk dot product and IncrementalPCA

   :param a: dask.Array of dataset a
   :param b: dask.Array of dataset b
   :param max_chunk_size: Chunk size for Incremental fit and transform, the larger the better as long as MEM is enough
   :param random_state:

   :returns:
   :rtype: Top CCA components


.. py:function:: calculate_direct_confusion(left_part, right_part)

   Given 2 dataframe for left/source and right/target dataset,
   calculate the direct confusion matrix based on co-cluster labels.
   Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
   The returned confusion matrix will be the form of source-cluster by target-cluster.

   :param left_part:
   :param right_part:


