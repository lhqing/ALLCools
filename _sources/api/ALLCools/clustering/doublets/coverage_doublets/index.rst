:py:mod:`ALLCools.clustering.doublets.coverage_doublets`
========================================================

.. py:module:: ALLCools.clustering.doublets.coverage_doublets

.. autoapi-nested-parse::

   This script implements the AMULET algorithm for doublet detection on snmC like data.

   The AMULET algorithm consider coverage information of each cell to predict doublets.
   It is complementary to the Scrublet like algorithm that uses molecular information.
   This algorithm is only suitable for DNA based sequencing data (and where input genome DNA copy is certain).

   Reference:
   Thibodeau, Asa, Alper Eroglu, Christopher S. McGinnis, Nathan Lawlor, Djamel Nehar-Belaid, Romy Kursawe,
   Radu Marches, et al. 2021. “AMULET: A Novel Read Count-Based Method for Effective Multiplet Detection from
   Single Nucleus ATAC-Seq Data.” Genome Biology 22 (1): 252.



Module Contents
---------------

.. py:function:: _calculate_cell_record(allc_path, output_path, cov_cutoff=2, resolution=100)

   Count the high coverage bins for each cell, save results to json


.. py:function:: calculate_blacklist_region(region_records, alpha=0.01)

   Collect highly covered regions by region-wise poisson FDR p value < alpha


.. py:function:: _calculate_cell_final_values(output_path, region_blacklist)

   Calculate final cell values while remove blacklist


.. py:function:: coverage_doublets(allc_dict: dict, resolution: int = 100, cov_cutoff=2, region_alpha=0.01, tmp_dir='doublets_temp_dir', cpu=1, keep_tmp=False)

   Quantify cell high coverage bins for doublets evaluation

   :param allc_dict: dict with cell_id as key, allc_path as value
   :param resolution: genome bin resolution to quantify, bps
   :param cov_cutoff: cutoff the cov, sites within cov_cutoff < cov <= 2 * cov_cutoff will be count
   :param region_alpha: FDR adjusted P-value cutoff
   :param tmp_dir: temporary dir to save the results
   :param cpu: number of cpu to use
   :param keep_tmp: Whether save the tem_dir for debugging


