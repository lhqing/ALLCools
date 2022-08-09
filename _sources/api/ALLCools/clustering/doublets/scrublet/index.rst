:py:mod:`ALLCools.clustering.doublets.scrublet`
===============================================

.. py:module:: ALLCools.clustering.doublets.scrublet


Module Contents
---------------

.. py:class:: MethylScrublet(sim_doublet_ratio=2.0, n_neighbors=None, expected_doublet_rate=0.1, stdev_doublet_rate=0.02, metric='euclidean', random_state=0, n_jobs=-1)

   Methods for calling doublets using the Scrublet algorithm on methylation dataset.

   .. py:method:: fit(mc, cov, clusters=None, batches=None)

      Fit the model to predict doublets.


   .. py:method:: simulate_doublets()

      Simulate doublets by adding the counts of random observed cell pairs.


   .. py:method:: pca()

      Perform PCA on data.


   .. py:method:: get_knn_graph(data)

      Get nearest neighbor graph for data.


   .. py:method:: calculate_doublet_scores()

      Calculate doublet scores for observed and simulated data.


   .. py:method:: call_doublets(threshold=None)

      Predicts doublets based on doublet scores.


   .. py:method:: plot()

      Plot the doublet score ROC curve and the doublet score histogram.


   .. py:method:: _plot_cluster_dist()



