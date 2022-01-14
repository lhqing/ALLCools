:py:mod:`ALLCools.clustering`
=============================

.. py:module:: ALLCools.clustering

.. autoapi-nested-parse::

   Basic Cellular Analysis Functions



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
   balanced_pca/index.rst
   dmg/index.rst
   incremental_pca/index.rst
   integration/index.rst
   mcad/index.rst
   pvclust/index.rst
   scrublet/index.rst


Package Contents
----------------

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


.. py:function:: tsne(adata, obsm='X_pca', metric: Union[str, Callable] = 'euclidean', exaggeration: float = -1, perplexity: int = 30, n_jobs: int = -1)

   Calculating T-SNE embedding with the openTSNE package :cite:p:`Policar2019` and
   parameter optimization strategy described in :cite:p:`Kobak2019`.

   :param adata: adata object with principle components or equivalent matrix stored in .obsm
   :param obsm: name of the matrix in .obsm that can be used as T-SNE input
   :param metric: Any metric allowed by PyNNDescent (default: 'euclidean')
   :param exaggeration: The exaggeration to use for the embedding
   :param perplexity: The perplexity to use for the embedding
   :param n_jobs: Number of CPUs to use

   :returns:
   :rtype: T-SNE embedding will be stored at adata.obsm["X_tsne"]


.. py:function:: balanced_pca(adata: anndata.AnnData, groups: str = 'pre_clusters', max_cell_prop=0.1, n_comps=200, scale=False)

   Given a categorical variable (e.g., a pre-clustering label), perform balanced PCA by downsample
   cells in the large categories to make the overall population more balanced, so the PCs are expected
   to represent more variance among small categories.

   :param adata: adata after preprocessing and feature selection steps
   :param groups: the name of the categorical variable in adata.obsm
   :param max_cell_prop: any single category with cells > `n_cell * max_cell_prop` will be downsampled to this number.
   :param n_comps: Number of components in PCA
   :param scale: whether to scale the input matrix before PCA

   :returns:
   :rtype: adata with PC information stored in obsm, varm and uns like the :func:`scanpy.tl.pca` do.


.. py:function:: significant_pc_test(adata: anndata.AnnData, p_cutoff=0.1, update=True, obsm='X_pca', downsample=50000)

   Perform two-sample Kolmogorov-Smirnov test for goodness of fit on two adjacent PCs,
   select top PCs based on the `p_cutoff`. Top PCs have significantly different distributions, while
   later PCs only capturing random noise will have larger p-values. An idea from :cite:p:`Zeisel2018`.

   :param adata: adata with PC matrix calculated and stored in adata.obsm
   :param p_cutoff: the p-value cutoff to select top PCs
   :param update: Whether modify adata.obsm and only keep significant PCs
   :param obsm: name of the PC matrix in adata.obsm
   :param downsample: If the dataset is too large, downsample the cells before testing.

   :returns: number of PCs selected
   :rtype: n_components


.. py:function:: log_scale(adata, method='standard', with_mean=False, with_std=True, max_value=10, scaler=None)

   Perform log transform and then scale the cell-by-feature matrix

   :param adata: adata with normalized, unscaled cell-by-feature matrix
   :param method: the type of scaler to use:
                  'standard' for :class:`sklearn.preprocessing.StandardScaler`;
                  'robust' for :class:`sklearn.preprocessing.RobustScaler`.
   :param with_mean: Whether scale with mean center
   :param with_std: Whether scale the std
   :param max_value: Whether clip large values after scale
   :param scaler: A fitted sklearn scaler, if provided, will only use it to transform the adata.

   :returns:
   :rtype: adata.X is scaled in place, the fitted scaler object will be return if the `scaler` parameter is None.


.. py:function:: get_pc_centers(adata: anndata.AnnData, group: str, outlier_label=None, obsm='X_pca')

   Get the cluster centroids of the PC matrix

   :param adata: adata with cluster labels in obs and PC in obsm
   :param group: the name of cluster labels in adata.obs
   :param outlier_label: if there are outlier labels in obs, will exclude outliers before calculating centroids.
   :param obsm: the key of PC matrix in obsm

   :returns: a dataframe for cluster centroids by PC
   :rtype: pc_center


.. py:class:: ReproduciblePCA(scaler, mc_type, adata=None, pca_obj=None, pc_loading=None, var_names=None, max_value=10)

   .. py:method:: mcds_to_adata(self, mcds)

      Get adata from MCDS with only selected features

      :param mcds:
      :type mcds: Input raw count MCDS object

      :returns: Adata with per-cell normalized mC fraction and selected features
      :rtype: adata


   .. py:method:: scale(self, adata)

      Perform {func}`log_scale <ALLCools.clustering.balanced_pca.log_scale>` with fitted scaler

      :param adata: Adata with per-cell normalized mC fraction and selected features

      :returns:
      :rtype: adata.X is transformed in place


   .. py:method:: pc_transform(self, adata)

      calculate the PC from adata.X and PC loading, store PCs in adata.obsm["X_pca"] and
      loadings in adata.varm["PCs"]

      :param adata: Adata with log_scale transformed mC fraction and selected features

      :returns:
      :rtype: PC information stored in adata.obsm and varm


   .. py:method:: mcds_to_adata_with_pc(self, mcds)

      From raw count MCDS to adata object with PCs using fitted scaler and PC loadings.
      Steps include select features, calculate per-cell normalized mC fractions,
      log_scale transform the data with fitted scaler, and finally add PC matrix.

      :param mcds: Raw count MCDS

      :returns: anndata object with per-cell normalized, log scale transformed matrix in .X and PCs in
                adata.obsm["X_pca"] and PC loadings in adata.varm["PCs"]. The scale and PC are done with fitted model.
      :rtype: adata


   .. py:method:: dump(self, path)

      Save the ReproduciblePCA to path



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

   .. py:method:: fit_predict(self, x, groups, obs_dim='cell', var_dim='gene', outlier='Outlier', cleanup=True, selected_pairs: List[tuple] = None)

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



.. py:function:: one_vs_rest_dmg(cell_meta, group, mcds=None, mcds_paths=None, obs_dim='cell', var_dim='gene', mc_type='CHN', top_n=1000, adj_p_cutoff=0.01, fc_cutoff=0.8, auroc_cutoff=0.8, max_cluster_cells=2000, max_other_fold=5, cpu=1)

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


.. py:function:: cluster_enriched_features(adata: anndata.AnnData, cluster_col: str, top_n=200, alpha=0.05, stat_plot=True, method='mc')

   Calculate top Cluster Enriched Features (CEF) from per-cell normalized dataset.
   An post-clustering feature selection step adapted from :cite:p:`Zeisel2018,La_Manno2021`
   and their great [cytograph2](https://github.com/linnarsson-lab/cytograph2) package.
   For details about CEF calculation, read the methods of :cite:p:`Zeisel2018`. Note that
   in original paper, they look for cluster-specific highly expressed genes as CEFs;
   for methylation, we are looking for hypo-methylation as CEFs, so the score and test is reversed.

   :param adata: adata containing per-cell normalized values.
                 For methylation fraction, the value need to be 1-centered
                 (1 means cell's average methylation), like those produced by
                 :func:`ALLCools.mcds.mcds.MCDS.add_mc_frac` with `normalize_per_cell=True`.
                 For RNA and ATAC, you can use per cell normalized counts.
                 Do not log transform the data before running this function
   :param cluster_col: The name of categorical variable in adata.obs
   :param top_n: Select top N CEFs for each cluster
   :param alpha: FDR corrected q-value cutoff
   :param stat_plot: Whether making some summary plots for the CEF calculation
   :param method: "mc" for methylation CEF (look for hypo-methylation),
                  "rna" and "atac" for the RNA and ATAC or any count based data
                  (use the original cytograph algorithm, look for higher value)

   :returns: * *Modify adata inplace, adding a dictionary in adata.uns called f"{cluster_col}_feature_enrichment"*
             * *The dictionary contains "qvals" (np.ndarray cluster-by-feature enrichment score q-value) and*
             * *"cluster_order" (cluster order of the "qvals")*


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


