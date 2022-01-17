:py:mod:`ALLCools.clustering.doublets`
======================================

.. py:module:: ALLCools.clustering.doublets


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   coverage_doublets/index.rst
   scrublet/index.rst


Package Contents
----------------

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



.. py:function:: coverage_doublets(allc_dict: dict, resolution: int = 100, cov_cutoff=2, region_alpha=0.01, tmp_dir='doublets_temp_dir', cpu=1, keep_tmp=False)

   Quantify cell high coverage bins for doublets evaluation

   :param allc_dict: dict with cell_id as key, allc_path as value
   :param resolution: genome bin resolution to quantify, bps
   :param cov_cutoff: cutoff the cov, sites within cov_cutoff < cov <= 2 * cov_cutoff will be count
   :param region_alpha: FDR adjusted P-value cutoff
   :param tmp_dir: temporary dir to save the results
   :param cpu: number of cpu to use
   :param keep_tmp: Whether save the tem_dir for debugging


