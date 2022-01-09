:py:mod:`ALLCools.pseudo_cell.pseudo_cell_kmeans`
=================================================

.. py:module:: ALLCools.pseudo_cell.pseudo_cell_kmeans


Module Contents
---------------

.. py:function:: _kmeans_division(matrix, cells, max_pseudo_size, max_k=50)


.. py:function:: _calculate_pseudo_group(clusters, total_matrix, cluster_size_cutoff=100, max_pseudo_size=25)


.. py:function:: _merge_pseudo_cell(adata, aggregate_func, pseudo_group_key)


.. py:function:: generate_pseudo_cells(adata, cluster_col='leiden', obsm='X_pca', cluster_size_cutoff=100, max_pseudo_size=25, aggregate_func='downsample')

   Balance the clusters by merge or downsample cells within each cluster.
   We first group the data by pre-defined clusters (cluster_col),
   then run k-means clustering iteratively on clusters with size > cluster_size_cutoff,
   the k-means clusters are called cell groups, and the maximum cell group size < max_pseudo_size,
   Finally, we generate a new adata for the balanced dataset.

   :param adata: Original AnnData object, raw count in X is recommended if aggregate_func is sum.
   :param cluster_col: The clustering label for downsample
   :param obsm: The obsm key name to use for performing k-means clustering within clusters.
   :param cluster_size_cutoff: Cluster size smaller than the cutoff will not be downsample or aggregated.
   :param max_pseudo_size: Maximum number of cells in one pseudo-cell group
   :param aggregate_func: 'downsample' means randomly select one cell from one pseudo-cell group;
                          'sum' means sum up all values in a pseudo-cell group
                          'mean' means take the average of each feature in a pseudo-cell group
                          'median' means take the median of each feature in a pseudo-cell group


