:py:mod:`ALLCools.integration.seurat_class`
===========================================

.. py:module:: ALLCools.integration.seurat_class


Module Contents
---------------

.. py:data:: CPU_COUNT
   

   

.. py:function:: top_features_idx(data, n_features)

   Select top features with the highest importance in CCs.

   :param data: data.shape = (n_cc, total_features)
   :param n_features: number of features to select

   :returns: **features_idx**
   :rtype: np.array


.. py:function:: find_neighbor(cc1, cc2, k, random_state=0, n_jobs=-1)

   Find all four way of neighbors for two datasets.

   :param cc1: cc for dataset 1
   :param cc2: cc for dataset 2
   :param k: number of neighbors
   :param random_state: random seed

   :rtype: 11, 12, 21, 22 neighbor matrix in shape (n_cell, k)


.. py:function:: find_mnn(G12, G21, kanchor)

   Calculate mutual nearest neighbor for two datasets.


.. py:function:: min_max(tmp, q_left=1, q_right=90)

   Normalize to q_left, q_right quantile to 0, 1, and cap extreme values.


.. py:function:: filter_anchor(anchor, adata_ref=None, adata_qry=None, high_dim_feature=None, k_filter=200, random_state=0, n_jobs=-1)

   Check if an anchor is still an anchor when only using the high_dim_features to construct KNN graph.

   If not, remove the anchor.


.. py:function:: score_anchor(anchor, G11, G12, G21, G22, k_score=30, Gp1=None, Gp2=None, k_local=50)

   Score the anchor by the number of shared neighbors.

   :param anchor: anchor in shape (n_anchor, 2)
   :param G11: neighbor graph of dataset 1
   :param G12: neighbor graph of dataset 1 to 2
   :param G21: neighbor graph of dataset 2 to 1
   :param G22: neighbor graph of dataset 2
   :param k_score: number of neighbors to score the anchor
   :param Gp1: Intra-dataset1 kNN graph
   :param Gp2: Intra-dataset2 kNN graph
   :param k_local: number of neighbors to calculate the local score

   :returns: **anchor with score in shape (n_anchor, 3)**
   :rtype: pd.DataFrame


.. py:function:: find_order(dist, ncell)

   Use dendrogram to find the order of dataset pairs.


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



