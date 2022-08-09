:py:mod:`ALLCools.clustering.ConsensusClustering`
=================================================

.. py:module:: ALLCools.clustering.ConsensusClustering


Module Contents
---------------

.. py:function:: _adata_to_coord_data(adata, coord_base)


.. py:function:: _r1_normalize(cmat)

   Perofrm R1 normalization on the confusion matrix.

   Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

   Normalize the confusion matrix based on the total number of cells in each class
   x(i,j) = max(cmat(i,j)/diagnol(i),cmat(j,i)/diagnol(j))
   confusion rate between i and j is defined by the maximum ratio i is confused as j or j is confused as i.

   :param cmat: the confusion matrix

   :rtype: the R1 normalized confusion matrix


.. py:function:: _r2_normalize(cmat)

   Perofrm R2 normalization on the confusion matrix.

   Adapted from https://github.com/SCCAF/sccaf/blob/develop/SCCAF/__init__.py

   Normalize the confusion matrix based on the total number of cells.
   x(i,j) = max(cmat(i,j)+cmat(j,i)/N)
   N is total number of cells analyzed.
   Confusion rate between i and j is defined by the sum of i confused as j or j confused as i.
   Then divide by total number of cells.

   :param cmat: the confusion matrix

   :rtype: the R2 normalized confusion matrix


.. py:function:: _leiden_runner(g, random_states, partition_type, **partition_kwargs)

   Run leiden on a graph and return the partition.

   The leiden clustering repeated len(random_states) times with different random states,
   return all clusters as a pd.DataFrame.


.. py:function:: _split_train_test_per_group(x, y, frac, max_train, random_state)

   Split train test for each cluster and make sure there are enough cells for train.


.. py:function:: single_supervise_evaluation(clf, x_train, y_train, x_test, y_test, r1_norm_step=0.05, r2_norm_step=0.05)

   Run supervise evaluation on confusion matrix.


.. py:class:: ConsensusClustering(model=None, n_neighbors=25, metric='euclidean', min_cluster_size=10, leiden_repeats=200, leiden_resolution=1, target_accuracy=0.95, consensus_rate=0.7, random_state=0, train_frac=0.5, train_max_n=500, max_iter=50, n_jobs=-1)

   .. py:method:: add_data(x)


   .. py:method:: fit_predict(x, leiden_kwds=None)


   .. py:method:: compute_neighbors()

      Calculate KNN graph


   .. py:method:: multi_leiden_clustering(partition_type=None, partition_kwargs=None, use_weights=True, n_iterations=-1)

      Run multiple leiden clustering with different random seeds and summarize the results.


   .. py:method:: _summarize_multi_leiden()

      Summarize the multi_leiden results.

      Generate a raw cluster version simply based on the hamming distance
      between cells and split cluster with cutoff (consensus_rate)


   .. py:method:: _create_model(n_estimators=1000)

      Init default model


   .. py:method:: supervise_learning()

      Perform supervised learning and cluster merge process


   .. py:method:: final_evaluation()

      Evaluate the final model


   .. py:method:: save(output_path)

      Save the model


   .. py:method:: plot_leiden_cases(coord_data, coord_base='umap', plot_size=3, dpi=300, plot_n_cases=4, s=3)

      Show some leiden runs with the biggest different as measured by ARI


   .. py:method:: plot_before_after(coord_data, coord_base='umap', plot_size=3, dpi=300)

      Plot the raw clusters from multi-leiden and final clusters after merge


   .. py:method:: plot_steps(coord_data, coord_base='umap', plot_size=3, dpi=300)

      Plot the supervised learning and merge steps


   .. py:method:: plot_merge_process(plot_size=3)

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


