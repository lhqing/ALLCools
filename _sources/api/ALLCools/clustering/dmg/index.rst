:py:mod:`ALLCools.clustering.dmg`
=================================

.. py:module:: ALLCools.clustering.dmg


Module Contents
---------------

.. py:function:: _single_pairwise_dmg(cluster_l, cluster_r, top_n, adj_p_cutoff, delta_rate_cutoff, auroc_cutoff, adata_dir, dmg_dir)

   Calculate DMG between a pair of adata file


.. py:class:: PairwiseDMG(max_cell_per_group=1000, top_n=10000, adj_p_cutoff=0.001, delta_rate_cutoff=0.3, auroc_cutoff=0.9, random_state=0, n_jobs=1, verbose=True)

   .. py:method:: fit_predict(self, x, groups, var_dim, obs_dim='cell', outlier='Outlier', cleanup=True, selected_pairs: List[tuple] = None)

      provide data and perform the pairwise DMG

      :param x: 2D cell-by-feature xarray.DataArray
      :param groups: cluster labels
      :param obs_dim: name of the cell dim
      :param var_dim: name of the feature dim
      :param outlier: name of the outlier group, if provided, will ignore this label
      :param cleanup: Whether to delete the group adata file
      :param selected_pairs: By default, pairwise DMG will calculate all possible pairs between all the groups, which might be very
                             time consuming if the group number is large. With this parameter, you may provide a list of cluster pairs


   .. py:method:: _save_cluster_adata(self)

      Save each group into separate adata, this way reduce the memory during parallel


   .. py:method:: _pairwise_dmg(self)

      pairwise DMG runner, result save to self.dmg_table


   .. py:method:: _cleanup(self)

      Delete group adata files


   .. py:method:: aggregate_pairwise_dmg(self, adata, groupby, obsm='X_pca')

      Aggregate pairwise DMG results for each cluster, rank DMG for the cluster by the sum of
      AUROC * cluster_pair_similarity
      This way, the DMGs having large AUROC between similar clusters get more weights

      :param adata:
      :param groupby:
      :param obsm:



.. py:function:: _single_ovr_dmg(cell_label, mcds, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff, auroc_cutoff)

   single one vs rest DMG runner


.. py:function:: _one_vs_rest_dmr_runner(cell_meta, group, cluster, max_cluster_cells, max_other_fold, mcds_paths, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff, auroc_cutoff, verbose=True)

   one vs rest DMG runner


.. py:function:: one_vs_rest_dmg(cell_meta, group, mcds=None, mcds_paths=None, obs_dim='cell', var_dim='gene', mc_type='CHN', top_n=1000, adj_p_cutoff=0.01, fc_cutoff=0.8, auroc_cutoff=0.8, max_cluster_cells=2000, max_other_fold=5, cpu=1, verbose=True)

   Calculating cluster marker genes using one-vs-rest strategy.

   :param cell_meta: cell metadata containing cluster labels
   :param group: the name of the cluster label column
   :param mcds: cell-by-gene MCDS object for calculating DMG. Provide either mcds_paths or mcds.
   :param mcds_paths: cell-by-gene MCDS paths for calculating DMG. Provide either mcds_paths or mcds.
   :param obs_dim: dimension name of the cells
   :param var_dim: dimension name of the features
   :param mc_type: value to select methylation type in the mc_type dimension
   :param top_n: report top N DMGs
   :param adj_p_cutoff: adjusted P value cutoff to report significant DMG
   :param fc_cutoff: mC fraction fold change cutoff to report significant DMG
   :param auroc_cutoff: AUROC cutoff to report significant DMG
   :param max_cluster_cells: The maximum number of cells from a group, downsample large group to this number
   :param max_other_fold: The fold of other cell numbers comparing
   :param cpu: number of cpus

   :returns: pandas Dataframe of the one-vs-rest DMGs
   :rtype: dmg_table


