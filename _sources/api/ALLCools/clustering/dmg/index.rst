:py:mod:`ALLCools.clustering.dmg`
=================================

.. py:module:: ALLCools.clustering.dmg


Module Contents
---------------

.. py:function:: single_pairwise_dmg(cluster_l, cluster_r, top_n, adj_p_cutoff, delta_rate_cutoff, auroc_cutoff, adata_dir, dmg_dir)


.. py:class:: PairwiseDMG(max_cell_per_group=1000, top_n=10000, adj_p_cutoff=0.001, delta_rate_cutoff=0.3, auroc_cutoff=0.9, random_state=0, n_jobs=1)

   .. py:method:: fit_predict(self, x, groups, obs_dim='cell', var_dim='gene', outlier='Outlier', cleanup=True, selected_pairs=None)


   .. py:method:: _save_cluster_adata(self)


   .. py:method:: _pairwise_dmg(self)


   .. py:method:: _cleanup(self)


   .. py:method:: aggregate_pairwise_dmg(self, adata, groupby, obsm='X_pca')



.. py:function:: single_ovr_dmg(cell_label, mcds, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff, auroc_cutoff)


.. py:function:: one_vs_rest_dmg(cell_meta, group, mcds=None, mcds_paths=None, obs_dim='cell', var_dim='gene', mc_type='CHN', top_n=1000, adj_p_cutoff=0.01, fc_cutoff=0.8, auroc_cutoff=0.8, max_cluster_cells=2000, max_other_fold=5, cpu=1)

   Calculating cluster marker genes using one-vs-rest strategy.

   :param cell_meta:
   :param group:
   :param mcds:
   :param mcds_paths:
   :param obs_dim:
   :param var_dim:
   :param mc_type:
   :param top_n:
   :param adj_p_cutoff:
   :param fc_cutoff:
   :param auroc_cutoff:
   :param max_cluster_cells:
   :param max_other_fold:


.. py:function:: _one_vs_rest_dmr_runner(cell_meta, group, cluster, max_cluster_cells, max_other_fold, mcds_paths, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff, auroc_cutoff)


