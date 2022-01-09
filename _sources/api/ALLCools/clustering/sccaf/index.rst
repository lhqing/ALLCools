:py:mod:`ALLCools.clustering.sccaf`
===================================

.. py:module:: ALLCools.clustering.sccaf


Module Contents
---------------

.. py:data:: color_long
   :annotation: = ['#e6194b', '#3cb44b', '#ffe119', '#0082c8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',...

   

.. py:function:: run_BayesianGaussianMixture(Y, K)

   For K-means clustering

   Input
   -----
   Y: the expression matrix
   K: number of clusters

   :returns:
   :rtype: clusters assigned to each cell.


.. py:function:: bhattacharyya_distance(repr1, repr2)

   Calculates Bhattacharyya distance (https://en.wikipedia.org/wiki/Bhattacharyya_distance).


.. py:function:: bhattacharyya_matrix(prob, flags=None)


.. py:function:: binary_accuracy(X, y, clf)


.. py:function:: normalize_confmat1(cmat, mod='1')

   Normalize the confusion matrix based on the total number of cells in each class
   x(i,j) = max(cmat(i,j)/diagnol(i),cmat(j,i)/diagnol(j))
   confusion rate between i and j is defined by the maximum ratio i is confused as j or j is confused as i.

   Input
   cmat: the confusion matrix

   :returns:
   :rtype: the normalized confusion matrix


.. py:function:: normalize_confmat2(cmat)

   Normalize the confusion matrix based on the total number of cells.
   x(i,j) = max(cmat(i,j)+cmat(j,i)/N)
   N is total number of cells analyzed.
   Confusion rate between i and j is defined by the sum of i confused as j or j confused as i.
   Then divide by total number of cells.

   Input
   cmat: the confusion matrix

   :returns:
   :rtype: the normalized confusion matrix


.. py:function:: cluster_adjmat(xmat, resolution=1, cutoff=0.1)

   Cluster the groups based on the adjacent matrix.
   Use the cutoff to discretize the matrix used to construct the adjacent graph.
   Then cluster the graph using the louvain clustering with a resolution value.
   As the adjacent matrix is binary, the default resolution value is set to 1.

   Input
   -----
   xmat: `numpy.array` or sparse matrix
       the reference matrix/normalized confusion matrix
   cutoff: `float` optional (default: 0.1)
       threshold used to binarize the reference matrix
   resolution: `float` optional (default: 1.0)
       resolution parameter for louvain clustering

   :returns:
   :rtype: new group names.


.. py:function:: msample(x, n, frac)

   sample the matrix by number or by fraction.
   if the fraction is larger than the sample number, use number for sampling. Otherwise, use fraction.

   Input
   -----
   x: the matrix to be split
   n: number of vectors to be sampled
   frac: fraction of the total matrix to be sampled

   :returns:
   :rtype: sampled selection.


.. py:function:: train_test_split_per_type(X, y, n=100, frac=0.8)

   This function is identical to train_test_split, but can split the data either based on number of cells or by fraction.

   Input
   -----
   X: `numpy.array` or sparse matrix
       the feature matrix
   y: `list of string/int`
       the class assignments
   n: `int` optional (default: 100)
       maximum number sampled in each label
   fraction: `float` optional (default: 0.8)
       Fraction of data included in the training set. 0.5 means use half of the data for training,
       if half of the data is fewer than maximum number of cells (n).

   :returns:
   :rtype: X_train, X_test, y_train, y_test


.. py:function:: SCCAF_assessment(*args, **kwargs)

   Assessment of clustering reliability using self-projection.
   It is the same as the self_projection function.


.. py:function:: self_projection(X, cell_types, classifier='LR', penalty='l1', sparsity=0.5, fraction=0.5, solver='liblinear', n=0, cv=5, whole=False, n_jobs=None)

   This is the core function for running self-projection.

   Input
   -----
   X: `numpy.array` or sparse matrix
       the expression matrix, e.g. ad.raw.X.
   cell_types: `list of String/int`
       the cell clustering assignment
   classifier: `String` optional (defatul: 'LR')
       a machine learning model in "LR" (logistic regression),         "RF" (Random Forest), "GNB"(Gaussion Naive Bayes), "SVM" (Support Vector Machine) and "DT"(Decision Tree).
   penalty: `String` optional (default: 'l2')
       the standardization mode of logistic regression. Use 'l1' or 'l2'.
   sparsity: `fload` optional (default: 0.5)
       The sparsity parameter (C in sklearn.linear_model.LogisticRegression) for the logistic regression model.
   fraction: `float` optional (default: 0.5)
       Fraction of data included in the training set. 0.5 means use half of the data for training,
       if half of the data is fewer than maximum number of cells (n).
   n: `int` optional (default: 100)
       Maximum number of cell included in the training set for each cluster of cells.
       only fraction is used to split the dataset if n is 0.
   cv: `int` optional (default: 5)
       fold for cross-validation on the training set.
       0 means no cross-validation.
   whole: `bool` optional (default: False)
       if measure the performance on the whole dataset (include training and test).
   n_jobs: `int` optional, number of threads to use with the different classifiers (default: None - unlimited).

   :returns: * *y_prob, y_pred, y_test, clf*
             * **y_prob** (`matrix of float`) -- prediction probability
             * **y_pred** (`list of string/int`) -- predicted clustering of the test set
             * **y_test** (`list of string/int`) -- real clustering of the test set
             * **clf** (*the classifier model.*)


.. py:function:: make_unique(dup_list)

   Make a name list unique by adding suffix "_%d". This function is identical to the make.unique function in R.

   Input
   -----
   dup_list: a list

   :returns:
   :rtype: a unique list with the same length as the input.


.. py:function:: confusion_matrix(y_test, y_pred, clf, labels=None)

   Get confusion matrix based on the test set.

   Input
   -----
   y_test, y_pred, clf: same as in self_projection

   :returns:
   :rtype: the confusion matrix


.. py:function:: per_cluster_accuracy(mtx, ad=None, clstr_name='louvain')

   Measure the accuracy of each cluster and put into a metadata slot.
   So the reliability of each cluster can be visualized.

   Input
   -----
   mtx: `pandas.dataframe`
       the confusion matrix
   ad: `AnnData`
       anndata object
   clstr_name: `String`
       the name of the clustering


.. py:function:: per_cell_accuracy(X, cell_types, clf)


.. py:function:: get_topmarkers(clf, names, topn=10)

   Get the top weighted features from the logistic regressioin model.

   Input
   -----
   clf: the logistic regression classifier
   names: `list of Strings`
       the names of the features (the gene names).
   topn: `int`
       number of top weighted featured to be returned.

   :returns:
   :rtype: list of markers for each of the cluster.


.. py:function:: eu_distance(X, gp1, gp2, cell)

   Measure the euclidean distance between two groups of cells and the third group.

   Input
   -----
   X: `np.array` or `sparse matrix`
       the total expression matrix
   gp1: `bool list`
       group1 of cells
   gp2: `bool list`
       group2 of cells
   cell: `bool list`
       group3 of cells, the group to be compared with gp1 and gp2.

   :returns: * `float value`
             * *the average distance difference.*


.. py:function:: get_distance_matrix(X, clusters, labels=None, metric='euclidean')

   Get the mean distance matrix between all clusters.

   Input
   -----
   X: `np.array` or `sparse matrix`
       the total expression matrix
   clusters: `string list`
       the assignment of the clusters
   labels: `string list`
       the unique labels of the clusters
   metric: `string` (optional, default: euclidean)
       distance metrics, see (http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html)

   :returns: the all-cluster to all-cluster distance matrix.
   :rtype: `np.array`


.. py:function:: merge_cluster(ad, old_id, new_id, groups)


.. py:function:: find_high_resolution(ad, resolution=4, n=100)


.. py:function:: get_connection_matrix(ad_obs, key1, key2)


.. py:function:: SCCAF_optimize_all(adata, start_groups=None, min_acc=0.9, r1_norm_cutoff=0.5, r2_norm_cutoff=0.05, R1norm_step=0.01, R2norm_step=0.001, min_iter=3, max_iter=10, *args, **kwargs)

   adata: `AnnData`
       The AnnData object of the expression profile.
   min_acc: `float` optional (default: 0.9)
       The minimum self-projection accuracy to be optimized for.
       e.g., 0.9 means the clustering optimization (merging process)
       will not stop until the self-projection accuracy is above 90%.
   R1norm_cutoff: `float` optional (default: 0.5)
       The start cutoff for R1norm of confusion matrix.
       e.g., 0.5 means use 0.5 as a cutoff to discretize the confusion matrix after R1norm.
       the discretized matrix is used to construct the connection graph for clustering optimization.
   R2norm_cutoff: `float` optional (default: 0.05)
       The start cutoff for R2norm of confusion matrix.
   R1norm_step: `float` optional (default: 0.01)
       The reduce step for minimum R1norm value.
       Each round of optimization calls the function `SCCAF_optimize`.
       The start of the next round of optimization is based on a new
       cutoff for the R1norm of the confusion matrix. This cutoff is
       determined by the minimum R1norm value in the previous round minus the R1norm_step value.
   R2norm_step: `float` optional (default: 0.001)
       The reduce step for minimum R2norm value.


.. py:function:: SCCAF_optimize(ad, prefix='L1', use='raw', use_projection=False, R1norm_only=False, R2norm_only=False, dist_only=False, dist_not=True, plot=True, basis='umap', plot_dist=False, plot_cmat=False, mod='1', low_res=None, c_iter=3, n_iter=10, n_jobs=None, start_iter=0, sparsity=0.5, n=100, fraction=0.5, r1_norm_cutoff=0.1, r2_norm_cutoff=1, dist_cutoff=8, classifier='LR', mplotlib_backend=None, min_acc=1)

   This is a self-projection confusion matrix directed cluster optimization function.

   Input
   -----
   ad: `AnnData`
       The AnnData object of the expression profile.
   prefix: `String`, optional (default: 'L1')
       The name of the optimization, which set as a prefix.
       e.g., the prefix = 'L1', the start round of optimization clustering is based on
       'L1_Round0'. So we need to assign an over-clustering state as a start point.
       e.g., ad.obs['L1_Round0'] = ad.obs['louvain']
   use: `String`, optional (default: 'raw')
       Use what features to train the classifier. Three choices:
       'raw' uses all the features;
       'hvg' uses the highly variable genes in the anndata object ad.var_names slot;
       'pca' uses the PCA data in the anndata object ad.obsm['X_pca'] slot.
   R1norm_only: `bool` optional (default: False)
       If only use the confusion matrix(R1norm) for clustering optimization.
   R2norm_only: `bool` optional (default: False)
       If only use the confusion matrix(R2norm) for clustering optimization.
   dist_only: `bool` optional (default: False)
       If only use the distance matrix for clustering optimization.
   dist_not: `bool` optional (default: True)
       If not use the distance matrix for clustering optimization.
   plot: `bool` optional (default: True)
       If plot the self-projectioin results, ROC curves and confusion matrices,
       during the optimization.
   plot_tsne: `bool` optional (default: False)
       If plot the self-projectioin results as tSNE. If False, the results are plotted as UMAP.
   plot_dist: `bool` optional (default: False)
       If make a scatter plot of the distance compared with the confusion rate for each of the cluster.
   plot_cmat: `bool` optional (default: False)
       plot the confusion matrix or not.
   mod: `string` optional (default: '1')
       two directions of normalization of confusion matrix for R1norm.
   c_iter: `int` optional (default: 3)
       Number of iterations of sampling for the confusion matrix.
       The minimum value of confusion rate in all the iterations is used as the confusion rate between two clusters.
   n_iter： `int` optional (default: 10)
       Maximum number of iterations(Rounds) for the clustering optimization.
   start_iter： `int` optional (default: 0)
       The start round of the optimization. e.g., start_iter = 3,
       the optimization will start from ad.obs['%s_3'%prefix].
   sparsity: `fload` optional (default: 0.5)
       The sparsity parameter (C in sklearn.linear_model.LogisticRegression) for the logistic regression model.
   n: `int` optional (default: 100)
       Maximum number of cell included in the training set for each cluster of cells.
   n_jobs: `int` number of jobs/threads to use (default: None - unlimited).
   fraction: `float` optional (default: 0.5)
       Fraction of data included in the training set. 0.5 means use half of the data for training,
       if half of the data is fewer than maximum number of cells (n).
   R1norm_cutoff: `float` optional (default: 0.1)
       The cutoff for the confusion rate (R1norm) between two clusters.
       0.1 means we allow maximum 10% of the one cluster confused as another cluster.
   R2norm_cutoff: `float` optional (default: 1.0)
       The cutoff for the confusion rate (R2norm) between two clusters.
       1.0 means the confusion between any two cluster should not exceed 1% of the total number of cells.
   dist_cutoff: `float` optional (default: 8.0)
       The cutoff for the euclidean distance between two clusters of cells.
       8.0 means the euclidean distance between two cell types should be greater than 8.0.
   low_res: `str` optional
               the clustering boundary for under-clustering. Set a low resolution in louvain/leiden clustering and give
               the key as the underclustering boundary.
   classifier: `String` optional (default: 'LR')
       a machine learning model in "LR" (logistic regression),         "RF" (Random Forest), "GNB"(Gaussion Naive Bayes), "SVM" (Support Vector Machine) and "DT"(Decision Tree).
   mplotlib_backend: `matplotlib.backends.backend_pdf` optional
       MatPlotLib multi-page backend object instance, previously initialised (currently the only type supported is
       PdfPages).
   min_acc: `float`
               the minimum total accuracy to be achieved. Above this threshold, the optimization will stop.

   :returns: assigned as the clustering optimization results.
   :rtype: The modified anndata object, with a slot "%s_result"%prefix


