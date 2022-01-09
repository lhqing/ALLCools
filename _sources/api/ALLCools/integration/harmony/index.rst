:py:mod:`ALLCools.integration.harmony`
======================================

.. py:module:: ALLCools.integration.harmony


Module Contents
---------------

.. py:data:: logger
   

   

.. py:data:: ch
   

   

.. py:data:: formatter
   

   

.. py:function:: leiden_centroids(pcs, n_pcs=None, n_neighbors=15, resolution=1.5)

   Instead of kmeans, using leiden algorithm to find centroids


.. py:function:: run_harmony(data_mat: numpy.ndarray, meta_data: pandas.DataFrame, vars_use, theta=None, lamb=None, sigma=0.1, nclust=None, tau=0, block_size=0.05, max_iter_harmony=10, max_iter_kmeans=20, epsilon_cluster=1e-05, epsilon_harmony=0.0001, verbose=True, random_state=0, init_method='kmeans', n_pcs=None, n_neighbors=15, resolution=1.5, leiden_input='origin')

   Run Harmony.



.. py:class:: Harmony(Z, Phi, Phi_moe, Pr_b, sigma, theta, max_iter_harmony, max_iter_kmeans, epsilon_kmeans, epsilon_harmony, K, block_size, lamb, verbose, init_method, n_pcs, n_neighbors, resolution, leiden_input)

   Bases: :py:obj:`object`

   .. py:method:: result(self)


   .. py:method:: allocate_buffers(self)


   .. py:method:: init_cluster(self)


   .. py:method:: compute_objective(self)


   .. py:method:: harmonize(self, iter_harmony=10, verbose=True)


   .. py:method:: cluster(self)


   .. py:method:: update_R(self)


   .. py:method:: check_convergence(self, i_type)



.. py:function:: safe_entropy(x: numpy.array)


.. py:function:: moe_correct_ridge(Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb)


