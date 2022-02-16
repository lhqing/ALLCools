:py:mod:`ALLCools.integration.cca`
==================================

.. py:module:: ALLCools.integration.cca


Module Contents
---------------

.. py:function:: simple_cca(adata, group_col, n_components=50, random_state=0)


.. py:function:: incremental_cca(a, b, max_chunk_size=10000, random_state=0)

   Perform Incremental CCA by chunk dot product and IncrementalPCA

   :param a: dask.Array of dataset a
   :param b: dask.Array of dataset b
   :param max_chunk_size: Chunk size for Incremental fit and transform, the larger the better as long as MEM is enough
   :param random_state:

   :returns:
   :rtype: Top CCA components


