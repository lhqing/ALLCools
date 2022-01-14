:py:mod:`ALLCools.clustering.ConsensusClustering`
=================================================

.. py:module:: ALLCools.clustering.ConsensusClustering


Module Contents
---------------

.. py:function:: _r1_normalize(cmat)

   Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

   Normalize the confusion matrix based on the total number of cells in each class
   x(i,j) = max(cmat(i,j)/diagnol(i),cmat(j,i)/diagnol(j))
   confusion rate between i and j is defined by the maximum ratio i is confused as j or j is confused as i.

   Input
   cmat: the confusion matrix

   :returns:
   :rtype: the normalized confusion matrix


.. py:function:: _r2_normalize(cmat)

   Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

   Normalize the confusion matrix based on the total number of cells.
   x(i,j) = max(cmat(i,j)+cmat(j,i)/N)
   N is total number of cells analyzed.
   Confusion rate between i and j is defined by the sum of i confused as j or j confused as i.
   Then divide by total number of cells.

   Input
   cmat: the confusion matrix

   :returns:
   :rtype: the normalized confusion matrix


.. py:function:: _leiden_runner(g, random_states, partition_type, **partition_kwargs)

   run leiden clustering len(random_states) times with different random states,
   return all clusters as a pd.DataFrame


.. py:function:: _split_train_test_per_group(x, y, frac, max_train, random_state)

   Split train test for each cluster and make sure there are enough cells for train


.. py:function:: single_supervise_evaluation(clf, x_train, y_train, x_test, y_test, r1_norm_step=0.05, r2_norm_step=0.05)

   A single fit and merge cluster step


.. py:class:: ConsensusClustering(model=None, n_neighbors=25, metric='euclidean', min_cluster_size=10, leiden_repeats=200, leiden_resolution=1, target_accuracy=0.95, consensus_rate=0.7, random_state=0, train_frac=0.5, train_max_n=500, max_iter=50, n_jobs=-1)

   .. py:method:: add_data(self, x)


   .. py:method:: fit_predict(self, x, leiden_kwds=None)


   .. py:method:: compute_neighbors(self)

      Calculate KNN graph


   .. py:method:: multi_leiden_clustering(self, partition_type=None, partition_kwargs=None, use_weights=True, n_iterations=-1)

      Modified from scanpy, perform Leiden clustering multiple times with different random states


   .. py:method:: _summarize_multi_leiden(self)

      Summarize the multi_leiden results,
      generate a raw cluster version simply based on the hamming distance
      between cells and split cluster with cutoff (consensus_rate)


   .. py:method:: _create_model(self, n_estimators=1000)

      Init default model


   .. py:method:: supervise_learning(self)

      Perform supervised learning and cluster merge process


   .. py:method:: final_evaluation(self)

      Final evaluation of the model and assign outliers


   .. py:method:: save(self, output_path)

      Save the model


   .. py:method:: plot_leiden_cases(self, coord_data, coord_base='umap', plot_size=3, dpi=300, plot_n_cases=4, s=3)

      Show some leiden runs with biggest different as measured by ARI


   .. py:method:: plot_before_after(self, coord_data, coord_base='umap', plot_size=3, dpi=300)

      Plot the raw clusters from multi-leiden and final clusters after merge


   .. py:method:: plot_steps(self, coord_data, coord_base='umap', plot_size=3, dpi=300)

      Plot the supervised learning and merge steps


   .. py:method:: plot_merge_process(self, plot_size=3)

      Plot the change of accuracy during merge



.. py:function:: select_confusion_pairs(true_label, predicted_label, ratio_cutoff=0.001)

   Select cluster pairs that are confusing (ratio_cutoff) between true and predicted labels

   :param true_label:
   :type true_label: true cell labels
   :param predicted_label:
   :type predicted_label: predicted cell labels
   :param ratio_cutoff:
   :type ratio_cutoff: ratio of clusters cutoff to define confusion

   :returns: list of cluster pair tuples
   :rtype: confused_pairs


