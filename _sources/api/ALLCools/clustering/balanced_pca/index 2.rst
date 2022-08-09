:py:mod:`ALLCools.clustering.balanced_pca`
==========================================

.. py:module:: ALLCools.clustering.balanced_pca

.. autoapi-nested-parse::

   Perform Balanced PCA as well as
   Reproducible PCA (A class allow user provide fitted model and just transform MCDS to adata with PCs)



Module Contents
---------------

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



