:py:mod:`ALLCools.clustering.doublets.scrublet`
===============================================

.. py:module:: ALLCools.clustering.doublets.scrublet


Module Contents
---------------

.. py:class:: MethylScrublet(sim_doublet_ratio=2.0, n_neighbors=None, expected_doublet_rate=0.1, stdev_doublet_rate=0.02, metric='euclidean', random_state=0, n_jobs=-1)

   .. py:method:: fit(self, mc, cov, clusters=None, batches=None)


   .. py:method:: simulate_doublets(self)

      Simulate doublets by adding the counts of random observed cell pairs.


   .. py:method:: pca(self)


   .. py:method:: get_knn_graph(self, data)


   .. py:method:: calculate_doublet_scores(self)


   .. py:method:: call_doublets(self, threshold=None)


   .. py:method:: plot(self)


   .. py:method:: _plot_cluster_dist(self)



