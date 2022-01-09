:py:mod:`ALLCools.clustering`
=============================

.. py:module:: ALLCools.clustering


Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   feature_selection/index.rst


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   ConsensusClustering/index.rst
   art_of_tsne/index.rst
   balanced_knn/index.rst
   balanced_pca/index.rst
   dmg/index.rst
   incremental_pca/index.rst
   integration/index.rst
   mcad/index.rst
   pvclust/index.rst
   scanorama/index.rst
   sccaf/index.rst
   scrublet/index.rst


Package Contents
----------------

.. py:class:: ConsensusClustering(model=None, n_neighbors=25, metric='euclidean', min_cluster_size=10, leiden_repeats=200, leiden_resolution=1, target_accuracy=0.95, consensus_rate=0.7, random_state=0, train_frac=0.5, train_max_n=500, max_iter=50, n_jobs=-1)

   .. py:method:: add_data(self, x)


   .. py:method:: fit_predict(self, x, leiden_kwds=None)


   .. py:method:: compute_neighbors(self)


   .. py:method:: multi_leiden_clustering(self, partition_type=None, partition_kwargs=None, use_weights=True, n_iterations=-1)

      Modified from scanpy


   .. py:method:: _summarize_multi_leiden(self)


   .. py:method:: _create_model(self, n_estimators=1000)


   .. py:method:: supervise_learning(self)


   .. py:method:: final_evaluation(self)


   .. py:method:: save(self, output_path)


   .. py:method:: plot_leiden_cases(self, coord_data, coord_base='umap', plot_size=3, dpi=300, plot_n_cases=4, s=3)


   .. py:method:: plot_before_after(self, coord_data, coord_base='umap', plot_size=3, dpi=300)


   .. py:method:: plot_steps(self, coord_data, coord_base='umap', plot_size=3, dpi=300)


   .. py:method:: plot_merge_process(self, plot_size=3)



.. py:function:: select_confusion_pairs(true_label, predicted_label, ratio_cutoff=0.001)


.. py:function:: tsne(adata, obsm='X_pca', metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1)


.. py:function:: balanced_pca(adata, groups='pre_clusters', max_cell_prop=0.1, n_comps=200, scale=False)


.. py:function:: significant_pc_test(adata, p_cutoff=0.1, update=True, obsm='X_pca', downsample=50000)

   :param adata:
   :param p_cutoff:
   :param update:
   :param obsm:
   :param downsample:


.. py:function:: log_scale(adata, method='standard', with_mean=False, with_std=True, max_value=10, scaler=None)


.. py:function:: get_pc_centers(adata, group, outlier_label=None, obsm='X_pca')


.. py:class:: ReproduciblePCA(scaler, mc_type, adata=None, pca_obj=None, pc_loading=None, var_names=None, max_value=10)

   .. py:method:: mcds_to_adata(self, mcds)


   .. py:method:: scale(self, adata)


   .. py:method:: pc_transform(self, adata)


   .. py:method:: mcds_to_adata_with_pc(self, mcds)


   .. py:method:: dump(self, path)



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



.. py:class:: Dendrogram(nboot=1000, method_dist='correlation', method_hclust='average', n_jobs=-1)

   .. py:method:: fit(self, data)

      :param data: The data is in obs-by-var form, row is obs.


   .. py:method:: save(self, output_path)



.. py:class:: PairwiseDMG(max_cell_per_group=1000, top_n=10000, adj_p_cutoff=0.001, delta_rate_cutoff=0.3, auroc_cutoff=0.9, random_state=0, n_jobs=1)

   .. py:method:: fit_predict(self, x, groups, obs_dim='cell', var_dim='gene', outlier='Outlier', cleanup=True, selected_pairs=None)


   .. py:method:: _save_cluster_adata(self)


   .. py:method:: _pairwise_dmg(self)


   .. py:method:: _cleanup(self)


   .. py:method:: aggregate_pairwise_dmg(self, adata, groupby, obsm='X_pca')



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


.. py:function:: cluster_enriched_features(adata, cluster_col, top_n=200, alpha=0.05, stat_plot=True, method='mc')

   :param adata:
   :param cluster_col:
   :param top_n:
   :param alpha:
   :param stat_plot:
   :param method: "mc" for methylation fraction,
                  "rna" and "atac" for the original algorithm for the RNA


.. py:function:: filter_regions(adata, hypo_cutoff=None)

   Filter regions based on their

   :param adata:
   :param hypo_cutoff: min number of cells that are hypo-methylated (1) in this region.
                       If None, will use adata.shape[0] * 0.003


.. py:function:: lsi(adata, scale_factor=100000, n_components=100, algorithm='arpack', obsm='X_pca', random_state=0, fit_size=None)

   Run TF-IDF on the binarized adata.X, followed by TruncatedSVD and then scale the components by svd.singular_values_

   :param adata:
   :param scale_factor:
   :param n_components:
   :param algorithm:
   :param obsm:
   :param random_state:
   :param fit_size: Ratio or absolute int value, use to downsample when fitting the SVD to speed up run time.


.. py:function:: binarize_matrix(adata, cutoff=0.95)

   Binarize adata.X with adata.X > cutoff

   :param adata: AnnData object whose X is survival function matrix
   :param cutoff: Cutoff to binarize the survival function

   :returns:
   :rtype: None


.. py:function:: remove_black_list_region(adata, black_list_path, f=0.2)

   Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

   :param adata:
   :param black_list_path: Path to the black list bed file
   :param f: Fraction of overlap when calling bedtools intersect

   :returns:
   :rtype: None


.. py:function:: calculate_direct_confusion(left_part, right_part)

   Given 2 dataframe for left/source and right/target dataset,
   calculate the direct confusion matrix based on co-cluster labels.
   Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
   The returned confusion matrix will be the form of source-cluster by target-cluster.

   :param left_part:
   :param right_part:


