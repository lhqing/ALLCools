:py:mod:`ALLCools.clustering.incremental_pca`
=============================================

.. py:module:: ALLCools.clustering.incremental_pca


Module Contents
---------------

.. py:function:: _normalize_per_cell(matrix, cell_sum)


.. py:class:: IncrementalPCA(n_components=100, sparse=False, normalize_per_cell=True, log1p=True, scale=True, **kwargs)

   .. py:method:: fit(self, ds, use_cells=None, use_features=None, chunk=500000, cell_sum=None, var_dim='gene', obs_dim='cell', load_chunk=None, random_shuffle=True)


   .. py:method:: transform(self, ds, use_cells=None, chunk=100000)


   .. py:method:: fit_transform(self, ds, use_cells=None, use_features=None, chunk=500000, cell_sum=None, var_dim='gene', obs_dim='cell', load_chunk=None, random_shuffle=True)



