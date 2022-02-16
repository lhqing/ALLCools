:py:mod:`ALLCools.integration.confusion`
========================================

.. py:module:: ALLCools.integration.confusion


Module Contents
---------------

.. py:function:: calculate_direct_confusion(left_part, right_part)

   Given 2 dataframe for left/source and right/target dataset,
   calculate the direct confusion matrix based on co-cluster labels.
   Each dataframe only contain 2 columns, first is original cluster, second is co-cluster.
   The returned confusion matrix will be the form of source-cluster by target-cluster.

   :param left_part:
   :param right_part:


