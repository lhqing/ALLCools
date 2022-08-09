:py:mod:`ALLCools.integration.confusion`
========================================

.. py:module:: ALLCools.integration.confusion


Module Contents
---------------

.. py:function:: calculate_direct_confusion(*args, **kwargs)


.. py:function:: _get_overlap_score(left_values, right_values)


.. py:function:: calculate_overlap_score(left_part, right_part)

   Calculate the overlap score between intra-dataset clusters using co-cluster information.

   Input are 2 dataframes for left/source and right/target dataset,
   Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
   The returned confusion matrix will be the form of source-cluster by target-cluster.

   :param left_part: Dataframe for left/source dataset.
   :param right_part: Dataframe for right/target dataset.


.. py:function:: calculate_diagonal_score(confusion_matrix, col_group, row_group)

   Given a confusion matrix, evaluate the overall integration performance with the diagonal score.

   :param confusion_matrix: A confusion matrix.
   :param col_group: Integration group for the columns.
   :param row_group: Integration group for the rows.

   :rtype: float


.. py:function:: confusion_matrix_clustering(confusion_matrix, min_value=0, max_value=0.9, seed=0)

   Given a confusion matrix, bi-clustering the matrix using Leiden Algorithm.

   :param confusion_matrix: A confusion matrix. Row is query, column is reference.
   :param min_value: minimum value to be used as an edge weight.
   :param max_value: maximum value to be used as an edge weight. Larger value will be capped to this value.
   :param seed: random seed for Leiden Algorithm.

   :rtype: query_group, ref_group, ordered confusion_matrix


