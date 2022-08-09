:py:mod:`ALLCools.clustering.lsi`
=================================

.. py:module:: ALLCools.clustering.lsi


Module Contents
---------------

.. py:function:: tf_idf(data, scale_factor=100000, idf=None)


.. py:function:: lsi(adata, scale_factor=100000, n_components=100, algorithm='arpack', obsm='X_pca', random_state=0, fit_size=None)

   Run TF-IDF on the binarized adata.X, followed by TruncatedSVD and then scale the components by singular values.

   :param adata: AnnData object
   :param scale_factor: scale factor for TF-IDF
   :param n_components: number of components to keep
   :param algorithm: algorithm to use for TruncatedSVD
   :param obsm: key in adata.obsm to store the components in
   :param random_state: random state for reproducibility
   :param fit_size: Ratio or absolute int value, use to downsample when fitting the SVD to speed up run time.


.. py:class:: LSI(scale_factor=100000, n_components=100, algorithm='arpack', random_state=0, idf=None, model=None)

   .. py:method:: _downsample_data(data, downsample)


   .. py:method:: _get_data(data)
      :staticmethod:


   .. py:method:: fit(data, downsample=None)


   .. py:method:: fit_transform(data, downsample=None, obsm_name='X_lsi')


   .. py:method:: transform(data, chunk_size=50000, obsm_name='X_lsi')


   .. py:method:: save(path)



