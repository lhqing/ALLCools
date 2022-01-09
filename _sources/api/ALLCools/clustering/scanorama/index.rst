:py:mod:`ALLCools.clustering.scanorama`
=======================================

.. py:module:: ALLCools.clustering.scanorama

.. autoapi-nested-parse::

   Allow to chose different metric

   Modified from Scanorama, see LICENSE
   https://github.com/brianhie/scanorama/blob/master/LICENSE



Module Contents
---------------

.. py:data:: ALPHA
   :annotation: = 0.1

   

.. py:data:: APPROX
   :annotation: = True

   

.. py:data:: BATCH_SIZE
   :annotation: = 5000

   

.. py:data:: DIMRED
   :annotation: = 100

   

.. py:data:: HVG
   

   

.. py:data:: KNN
   :annotation: = 20

   

.. py:data:: N_ITER
   :annotation: = 500

   

.. py:data:: PERPLEXITY
   :annotation: = 1200

   

.. py:data:: SIGMA
   :annotation: = 15

   

.. py:data:: VERBOSE
   :annotation: = 2

   

.. py:function:: batch_correct_pc(adata, batch_series, correct=False, n_components=30, sigma=25, alpha=0.1, knn=30, metric='angular', **scanorama_kws)

   Batch correction PCA based on integration

   :param adata: one major adata
   :param batch_series: batch_series used for splitting adata
   :param correct: if True, adata.X will be corrected inplace, otherwise only corrected PCs are added to adata.obsm['X_pca']
   :param n_components: number of components in PCA
   :param sigma: Correction smoothing parameter on Gaussian kernel.
   :param alpha: Alignment score minimum cutoff.
   :param knn: Number of nearest neighbors to use for matching.
   :param metric: Metric to use in calculating KNN
   :param scanorama_kws: Other Parameters passed to integration function

   :returns:
   :rtype: adata


.. py:function:: correct(datasets_full, genes_list, return_dimred=False, batch_size=BATCH_SIZE, verbose=VERBOSE, ds_names=None, dimred=DIMRED, approx=APPROX, sigma=SIGMA, alpha=ALPHA, knn=KNN, return_dense=False, hvg=None, union=False, geosketch=False, geosketch_max=20000, seed=0, metric='manhattan')

   Integrate and batch correct a list of data sets.

   :param datasets_full: Data sets to integrate and correct.
   :type datasets_full: `list` of `scipy.sparse.csr_matrix` or of `numpy.ndarray`
   :param genes_list: List of genes for each data set.
   :type genes_list: `list` of `list` of `string`
   :param return_dimred: In addition to returning batch corrected matrices, also returns
                         integrated low-dimesional embeddings.
   :type return_dimred: `bool`, optional (default: `False`)
   :param batch_size: The batch size used in the alignment vector computation. Useful when
                      correcting very large (>100k samples) data sets. Set to large value
                      that runs within available memory.
   :type batch_size: `int`, optional (default: `5000`)
   :param verbose: When `True` or not equal to 0, prints logging output.
   :type verbose: `bool` or `int`, optional (default: 2)
   :param ds_names: When `verbose=True`, reports data set names in logging output.
   :type ds_names: `list` of `string`, optional
   :param dimred: Dimensionality of integrated embedding.
   :type dimred: `int`, optional (default: 100)
   :param approx: Use approximate nearest neighbors, greatly speeds up matching runtime.
   :type approx: `bool`, optional (default: `True`)
   :param sigma: Correction smoothing parameter on Gaussian kernel.
   :type sigma: `float`, optional (default: 15)
   :param alpha: Alignment score minimum cutoff.
   :type alpha: `float`, optional (default: 0.10)
   :param knn: Number of nearest neighbors to use for matching.
   :type knn: `int`, optional (default: 20)
   :param return_dense: Return `numpy.ndarray` matrices instead of `scipy.sparse.csr_matrix`.
   :type return_dense: `bool`, optional (default: `False`)
   :param hvg: Use this number of top highly variable genes based on dispersion.
   :type hvg: `int`, optional (default: None)
   :param seed: Random seed to use.
   :type seed: `int`, optional (default: 0)

   :returns: * *corrected, genes* -- By default (`return_dimred=False`), returns a two-tuple containing a
               list of `scipy.sparse.csr_matrix` each with batch corrected values,
               and a single list of genes containing the intersection of inputted
               genes.
             * *integrated, corrected, genes* -- When `return_dimred=False`, returns a three-tuple containing a list
               of `numpy.ndarray` with integrated low dimensional embeddings, a list
               of `scipy.sparse.csr_matrix` each with batch corrected values, and a
               a single list of genes containing the intersection of inputted genes.


.. py:function:: integrate(datasets_full, genes_list, batch_size=BATCH_SIZE, verbose=VERBOSE, ds_names=None, dimred=DIMRED, approx=APPROX, sigma=SIGMA, alpha=ALPHA, knn=KNN, geosketch=False, geosketch_max=20000, n_iter=1, union=False, hvg=None, seed=0, metric='manhattan')

   Integrate a list of data sets.

   :param datasets_full: Data sets to integrate and correct.
   :type datasets_full: `list` of `scipy.sparse.csr_matrix` or of `numpy.ndarray`
   :param genes_list: List of genes for each data set.
   :type genes_list: `list` of `list` of `string`
   :param batch_size: The batch size used in the alignment vector computation. Useful when
                      correcting very large (>100k samples) data sets. Set to large value
                      that runs within available memory.
   :type batch_size: `int`, optional (default: `5000`)
   :param verbose: When `True` or not equal to 0, prints logging output.
   :type verbose: `bool` or `int`, optional (default: 2)
   :param ds_names: When `verbose=True`, reports data set names in logging output.
   :type ds_names: `list` of `string`, optional
   :param dimred: Dimensionality of integrated embedding.
   :type dimred: `int`, optional (default: 100)
   :param approx: Use approximate nearest neighbors, greatly speeds up matching runtime.
   :type approx: `bool`, optional (default: `True`)
   :param sigma: Correction smoothing parameter on Gaussian kernel.
   :type sigma: `float`, optional (default: 15)
   :param alpha: Alignment score minimum cutoff.
   :type alpha: `float`, optional (default: 0.10)
   :param knn: Number of nearest neighbors to use for matching.
   :type knn: `int`, optional (default: 20)
   :param hvg: Use this number of top highly variable genes based on dispersion.
   :type hvg: `int`, optional (default: None)
   :param seed: Random seed to use.
   :type seed: `int`, optional (default: 0)

   :returns: Returns a two-tuple containing a list of `numpy.ndarray` with
             integrated low dimensional embeddings and a single list of genes
             containing the intersection of inputted genes.
   :rtype: integrated, genes


.. py:function:: correct_scanpy(adatas, **kwargs)

   Batch correct a list of `scanpy.api.AnnData`.

   :param adatas: Data sets to integrate and/or correct.
   :type adatas: `list` of `scanpy.api.AnnData`
   :param kwargs: See documentation for the `correct()` method for a full list of
                  parameters to use for batch correction.
   :type kwargs: `dict`

   :returns: * *corrected* -- By default (`return_dimred=False`), returns a list of new
               `scanpy.api.AnnData`.
             * *integrated, corrected* -- When `return_dimred=True`, returns a two-tuple containing a list of
               `np.ndarray` with integrated low-dimensional embeddings and a list
               of new `scanpy.api.AnnData`.


.. py:function:: integrate_scanpy(adatas, **kwargs)

   Integrate a list of `scanpy.api.AnnData`.

   :param adatas: Data sets to integrate.
   :type adatas: `list` of `scanpy.api.AnnData`
   :param kwargs: See documentation for the `integrate()` method for a full list of
                  parameters to use for batch correction.
   :type kwargs: `dict`

   :returns: Returns a list of `np.ndarray` with integrated low-dimensional
             embeddings.
   :rtype: integrated


.. py:function:: merge_datasets(datasets, genes, ds_names=None, verbose=True, union=False)


.. py:function:: check_datasets(datasets_full)


.. py:function:: reduce_dimensionality(X, dim_red_k=100)


.. py:function:: dimensionality_reduce(datasets, dimred=DIMRED)


.. py:function:: dispersion(X)


.. py:function:: process_data(datasets, genes, hvg=HVG, dimred=DIMRED, verbose=True)


.. py:function:: nn(ds1, ds2, knn=KNN, metric_p=2)


.. py:function:: nn_approx(ds1, ds2, knn=KNN, metric='manhattan', n_trees=10)


.. py:function:: fill_table(table, i, curr_ds, datasets, base_ds=0, knn=KNN, approx=APPROX, metric='manhattan')


.. py:data:: gs_idxs
   

   

.. py:function:: find_alignments_table(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE, prenormalized=False, geosketch=False, geosketch_max=20000, metric='manhattan')


.. py:function:: find_alignments(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE, alpha=ALPHA, prenormalized=False, geosketch=False, geosketch_max=20000, metric='manhattan')


.. py:function:: connect(datasets, knn=KNN, approx=APPROX, alpha=ALPHA, verbose=VERBOSE, metric='manhattan')


.. py:function:: handle_zeros_in_scale(scale, copy=True)

   Makes sure that whenever scale is zero, we handle it correctly.
   This happens in most scalers when we have constant features.
   Adapted from sklearn.preprocessing.data


.. py:function:: batch_bias(curr_ds, match_ds, bias, batch_size=None, sigma=SIGMA)


.. py:function:: transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma=SIGMA, cn=False, batch_size=None)


.. py:function:: assemble(datasets, verbose=VERBOSE, knn=KNN, sigma=SIGMA, approx=APPROX, alpha=ALPHA, expr_datasets=None, ds_names=None, batch_size=None, geosketch=False, geosketch_max=20000, alignments=None, matches=None, metric='manhattan')


.. py:function:: interpret_alignments(datasets, expr_datasets, genes, verbose=VERBOSE, knn=KNN, approx=APPROX, alpha=ALPHA, n_permutations=None, metric='manhattan')


