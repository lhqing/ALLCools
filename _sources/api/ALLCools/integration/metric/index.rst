:py:mod:`ALLCools.integration.metric`
=====================================

.. py:module:: ALLCools.integration.metric


Module Contents
---------------

.. py:function:: _purity(label1, label2)

   Calculate purity score.


.. py:function:: calculate_purity(adata, label1, label2)

   Calculate purity score of two kinds of labels for the same set of obs.

   :param adata: anndata object
   :param label1: column name of the first label
   :param label2: column name of the second label

   :rtype: purity score


.. py:function:: calculate_adjust_rand_index(adata, label1, label2)

   Calculate adjusted rand index of two kinds of labels for the same set of obs.

   :param adata: anndata object
   :param label1: column name of the first label
   :param label2: column name of the second label

   :rtype: adjusted rand index


.. py:function:: _alignment_score(obsm_list, k=20, random_state=0, downsample_obs=None, n_jobs=-1)

   Calculate alignment score of multiple datasets.


.. py:function:: calculate_alignment_score(adata, dataset_col, obsm_key, downsample_obs=None, k=20, random_state=0, n_jobs=-1)

   Calculate alignment score of multiple datasets.

   :param adata: anndata object
   :param dataset_col: adata.obs column name of the dataset label
   :param obsm_key: adata.obsm key of the obsm matrix
   :param downsample_obs: downsample the obs of each dataset to this size
   :param k: number of neighbors to use
   :param random_state: random state for KNN
   :param n_jobs: number of jobs to use for KNN

   :rtype: overall alignment score, alignment score per dataset


.. py:function:: calculate_cluster_alignment_score(adata, cluster_col, dataset_col, obsm_key, downsample_obs=None, k=20, random_state=0, n_jobs=-1)

   Calculate alignment score of multiple datasets for each cluster.

   :param adata: anndata object
   :param cluster_col: adata.obs column name of the cluster label
   :param dataset_col: adata.obs column name of the dataset label
   :param obsm_key: adata.obsm key of the obsm matrix
   :param downsample_obs: downsample the obs of each dataset to this size
   :param k: number of neighbors to use
   :param random_state: random state for KNN
   :param n_jobs: number of jobs to use for KNN

   :returns: * *overall alignment score for each cluster,*
             * *alignment score per dataset for each cluster*


.. py:function:: calculate_kbet_accept_rate(adata, dataset_col, obsm_key, k=20, test_size=1000, downsample_obs=5000)

   Calculate the average observed accept rate of KBET test.

   :param adata: anndata object
   :param dataset_col: adata.obs column name of the dataset label
   :param obsm_key: adata.obsm key of the obsm matrix
   :param k: number of neighbors to use, pass to k0 parameter of kBET
   :param test_size: number of samples to use for test, pass to testSize parameter of kBET
   :param downsample_obs: downsample the obs of each dataset to this size

   :returns: * *average observed accept rate of kBET test,*
             * *complete kBET test result in a dictionary*


