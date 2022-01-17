:py:mod:`ALLCools.clustering.ClusterMerging`
============================================

.. py:module:: ALLCools.clustering.ClusterMerging


Module Contents
---------------

.. py:class:: ClusterMerge(merge_criterion, stop_criterion=None, stop_clusters=-1, n_cells=200, metric='euclidean', method='average', label_concat_str='::')

   .. py:method:: _construct_tree(self)


   .. py:method:: _traverse(node, call_back)
      :staticmethod:


   .. py:method:: _merge_pair(self, pair, concat_str='::')


   .. py:method:: fit_predict(self, data_for_tree, cell_to_type, gene_mcds)



.. py:class:: PairwiseDMGCriterion(max_cell_per_group=100, top_n_markers=5, adj_p_cutoff=0.001, delta_rate_cutoff=0.3, auroc_cutoff=0.85, use_modality='either', random_state=0, n_jobs=10, verbose=False)

   .. py:method:: predict(self, pair_labels, pair_cells, pair_cell_types, pair_mcds, da_name='gene_da_frac')



