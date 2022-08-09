:py:mod:`ALLCools.clustering.incremental_pca`
=============================================

.. py:module:: ALLCools.clustering.incremental_pca


Module Contents
---------------

.. py:function:: _normalize_per_cell(matrix, cell_sum)

   Normalize matrix row sum to to cell_sum.


.. py:class:: IncrementalPCA(n_components=100, sparse=False, normalize_per_cell=True, log1p=True, scale=True, **kwargs)

   Incremental PCA to fit and transform large datasets.

   .. py:method:: fit(ds, use_cells=None, use_features=None, chunk=500000, cell_sum=None, var_dim='gene', obs_dim='cell', load_chunk=None, random_shuffle=True)

      Fit the dataset to get PC Loadings.


   .. py:method:: transform(ds, use_cells=None, chunk=100000)

      Transform the dataset to PCA space.


   .. py:method:: fit_transform(ds, use_cells=None, use_features=None, chunk=500000, cell_sum=None, var_dim='gene', obs_dim='cell', load_chunk=None, random_shuffle=True)

      Fit and transform the dataset.



