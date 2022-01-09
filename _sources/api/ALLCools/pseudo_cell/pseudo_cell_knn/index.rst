:py:mod:`ALLCools.pseudo_cell.pseudo_cell_knn`
==============================================

.. py:module:: ALLCools.pseudo_cell.pseudo_cell_knn


Module Contents
---------------

.. py:class:: ContractedExamplerSampler(data, n_components=30, normalize=False)

   .. py:method:: _sample_pulp_dist(self, n_kernels, pulp_size)


   .. py:method:: _select_dist_thresh(self, pulp_size, n_tests=100, pulp_thicken_ratio=1.2, robust_quantile=0.9)


   .. py:method:: _sample_fruit(self, n_kernels, pulp_size, max_iters, dist_thresh=None, ovlp_tol=0.2, min_pulp_size=None, k=100)


   .. py:method:: sample_contracted_examplers(self, n_examplers, n_neighbors, min_n_neighbors=None, ovlp_tol=0, dist_thresh=None, max_iters=100)



.. py:function:: sample_pseudo_cells(cell_meta, cluster_col, coords, target_pseudo_size, min_pseudo_size=None, ignore_small_cluster=False, n_components=30, pseudo_ovlp=0)


.. py:function:: generate_pseudo_cells(adata, cluster_col='leiden', obsm='X_pca', target_pseudo_size=100, min_pseudo_size=None, ignore_small_cluster=False, n_components=None, aggregate_func='downsample', pseudo_ovlp=0)


