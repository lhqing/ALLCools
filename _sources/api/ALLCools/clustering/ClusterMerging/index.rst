:py:mod:`ALLCools.clustering.ClusterMerging`
============================================

.. py:module:: ALLCools.clustering.ClusterMerging


Module Contents
---------------

.. py:class:: ClusterMerge(merge_criterion, stop_criterion=None, stop_clusters=-1, n_cells=200, metric='euclidean', method='average', label_concat_str='::')

   Perform cluster merge based on the given merge criterion.

   .. py:method:: _construct_tree()


   .. py:method:: _traverse(node, call_back)
      :staticmethod:


   .. py:method:: _merge_pair(pair, concat_str='::')


   .. py:method:: fit_predict(data_for_tree, cell_to_type, gene_mcds)

      Fit the model and predict the cluster merge.



.. py:class:: PairwiseDMGCriterion(max_cell_per_group=100, top_n_markers=5, adj_p_cutoff=0.001, delta_rate_cutoff=0.3, auroc_cutoff=0.85, use_modality='either', random_state=0, n_jobs=10, verbose=False)

   Perform pairwise DMG analysis to determine whether two clusters should be merged.

   .. py:method:: predict(pair_labels, pair_cells, pair_cell_types, pair_mcds, da_name='gene_da_frac')

      Predict whether two clusters should be merged.



