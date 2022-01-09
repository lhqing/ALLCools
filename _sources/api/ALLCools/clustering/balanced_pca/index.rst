:py:mod:`ALLCools.clustering.balanced_pca`
==========================================

.. py:module:: ALLCools.clustering.balanced_pca


Module Contents
---------------

.. py:function:: log_scale(adata, method='standard', with_mean=False, with_std=True, max_value=10, scaler=None)


.. py:function:: significant_pc_test(adata, p_cutoff=0.1, update=True, obsm='X_pca', downsample=50000)

   :param adata:
   :param p_cutoff:
   :param update:
   :param obsm:
   :param downsample:


.. py:function:: balanced_pca(adata, groups='pre_clusters', max_cell_prop=0.1, n_comps=200, scale=False)


.. py:function:: get_pc_centers(adata, group, outlier_label=None, obsm='X_pca')


.. py:class:: ReproduciblePCA(scaler, mc_type, adata=None, pca_obj=None, pc_loading=None, var_names=None, max_value=10)

   .. py:method:: mcds_to_adata(self, mcds)


   .. py:method:: scale(self, adata)


   .. py:method:: pc_transform(self, adata)


   .. py:method:: mcds_to_adata_with_pc(self, mcds)


   .. py:method:: dump(self, path)



