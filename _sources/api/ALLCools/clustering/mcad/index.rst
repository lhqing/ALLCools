:py:mod:`ALLCools.clustering.mcad`
==================================

.. py:module:: ALLCools.clustering.mcad


Module Contents
---------------

.. py:function:: remove_black_list_region(adata, black_list_path, f=0.2)

   Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

   :param adata:
   :param black_list_path: Path to the black list bed file
   :param f: Fraction of overlap when calling bedtools intersect

   :returns:
   :rtype: None


.. py:function:: binarize_matrix(adata, cutoff=0.95)

   Binarize adata.X with adata.X > cutoff

   :param adata: AnnData object whose X is survival function matrix
   :param cutoff: Cutoff to binarize the survival function

   :returns:
   :rtype: None


.. py:function:: filter_regions(adata, hypo_cutoff=None)

   Filter regions based on their

   :param adata:
   :param hypo_cutoff: min number of cells that are hypo-methylated (1) in this region.
                       If None, will use adata.shape[0] * 0.003


.. py:function:: tf_idf(data, scale_factor)


.. py:function:: lsi(adata, scale_factor=100000, n_components=100, algorithm='arpack', obsm='X_pca', random_state=0, fit_size=None)

   Run TF-IDF on the binarized adata.X, followed by TruncatedSVD and then scale the components by svd.singular_values_

   :param adata:
   :param scale_factor:
   :param n_components:
   :param algorithm:
   :param obsm:
   :param random_state:
   :param fit_size: Ratio or absolute int value, use to downsample when fitting the SVD to speed up run time.


