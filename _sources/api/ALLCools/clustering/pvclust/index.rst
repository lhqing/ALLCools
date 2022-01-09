:py:mod:`ALLCools.clustering.pvclust`
=====================================

.. py:module:: ALLCools.clustering.pvclust


Module Contents
---------------

.. py:function:: install_r_package(name)


.. py:function:: _hclust_to_scipy_linkage(result, plot=True)

   Turn R hclust result obj into scipy linkage matrix format


.. py:class:: Dendrogram(nboot=1000, method_dist='correlation', method_hclust='average', n_jobs=-1)

   .. py:method:: fit(self, data)

      :param data: The data is in obs-by-var form, row is obs.


   .. py:method:: save(self, output_path)



