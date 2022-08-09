:py:mod:`ALLCools.integration`
==============================

.. py:module:: ALLCools.integration


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cca/index.rst
   confusion/index.rst
   metric/index.rst
   seurat_class/index.rst


Package Contents
----------------

.. py:function:: calculate_diagonal_score(confusion_matrix, col_group, row_group)

   Given a confusion matrix, evaluate the overall integration performance with the diagonal score.

   :param confusion_matrix: A confusion matrix.
   :param col_group: Integration group for the columns.
   :param row_group: Integration group for the rows.

   :rtype: float


.. py:function:: calculate_direct_confusion(*args, **kwargs)


.. py:function:: calculate_overlap_score(left_part, right_part)

   Calculate the overlap score between intra-dataset clusters using co-cluster information.

   Input are 2 dataframes for left/source and right/target dataset,
   Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
   The returned confusion matrix will be the form of source-cluster by target-cluster.

   :param left_part: Dataframe for left/source dataset.
   :param right_part: Dataframe for right/target dataset.


.. py:function:: confusion_matrix_clustering(confusion_matrix, min_value=0, max_value=0.9, seed=0)

   Given a confusion matrix, bi-clustering the matrix using Leiden Algorithm.

   :param confusion_matrix: A confusion matrix. Row is query, column is reference.
   :param min_value: minimum value to be used as an edge weight.
   :param max_value: maximum value to be used as an edge weight. Larger value will be capped to this value.
   :param seed: random seed for Leiden Algorithm.

   :rtype: query_group, ref_group, ordered confusion_matrix


.. py:class:: SeuratIntegration(n_jobs=-1, random_state=0)

   Main class for Seurat integration.

   .. py:method:: _calculate_local_knn()

      Calculate local kNN graph for each dataset.

      If klocal is provided, we calculate the local knn graph to
      evaluate whether the anchor preserves local structure within the dataset.
      One can use a different obsm with key_local to compute knn for each dataset.


   .. py:method:: _get_all_pairs()


   .. py:method:: _prepare_matrix(i, j, key_anchor)


   .. py:method:: _calculate_mutual_knn_and_raw_anchors(i, j, U, V, k, k_anchor)

      Calculate the mutual knn graph and raw anchors.

      The results are saved to self.mutual_knn and self.raw_anchor.


   .. py:method:: find_anchor(adata_list, adata_names=None, k_local=None, key_local='X_pca', key_anchor='X', dim_red='pca', svd_algorithm='randomized', scale1=False, scale2=False, k_filter=None, n_features=200, n_components=None, max_cc_cells=50000, k_anchor=5, k_score=30, alignments=None)

      Find anchors for each dataset pair.


   .. py:method:: find_nearest_anchor(data, data_qry, ref, qry, key_correct='X_pca', npc=30, kweight=100, sd=1, random_state=0)

      Find the nearest anchors for each cell in data.


   .. py:method:: transform(data, ref, qry, key_correct, npc=30, k_weight=100, sd=1, chunk_size=50000, random_state=0, row_normalize=True)

      Transform query data to reference space.


   .. py:method:: integrate(key_correct, row_normalize=True, n_components=30, k_weight=100, sd=1, alignments=None)

      Integrate datasets by transform data matrices from query to reference data using the MNN information.


   .. py:method:: label_transfer(ref, qry, categorical_key=None, continuous_key=None, key_dist='X_pca', kweight=100, npc=30, sd=1, chunk_size=50000, random_state=0)

      Transfer labels from query to reference space.


   .. py:method:: save(output_path, save_local_knn=False, save_raw_anchor=False, save_mutual_knn=False)

      Save the model and results to disk.


   .. py:method:: load(input_path)
      :classmethod:

      Load integrator from file.


   .. py:method:: save_transfer_results_to_adata(adata, transfer_results, new_label_suffix='_transfer')
      :classmethod:

      Save transfer results to adata.



