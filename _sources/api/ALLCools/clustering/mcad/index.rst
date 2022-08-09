:py:mod:`ALLCools.clustering.mcad`
==================================

.. py:module:: ALLCools.clustering.mcad


Module Contents
---------------

.. py:function:: remove_black_list_region(adata, black_list_path, f=0.2)

   Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path.

   :param adata: AnnData object
   :param black_list_path: Path to the black list bed file
   :param f: Fraction of overlap when calling bedtools intersect


.. py:function:: remove_chromosomes(adata, exclude_chromosomes=None, include_chromosomes=None, chrom_col='chrom')

   Remove chromosomes from adata.var.


.. py:function:: binarize_matrix(adata, cutoff=0.95)

   Binarize adata.X with adata.X > cutoff

   :param adata: AnnData object whose X is survival function matrix
   :param cutoff: Cutoff to binarize the survival function

   :rtype: None


.. py:function:: filter_regions(adata, hypo_percent=0.5, n_cell=None, zscore_abs_cutoff=None)

   Filter regions based on % of cells having non-zero scores.

   :param adata: AnnData object
   :param hypo_percent: min % of cells that are non-zero in this region. If n_cell is provided, this parameter will be ignored.
   :param n_cell: number of cells that are non-zero in this region.
   :param zscore_abs_cutoff: absolute feature non-zero cell count zscore cutoff to remove lowest and highest coverage features.


