:py:mod:`ALLCools.dmr`
======================

.. py:module:: ALLCools.dmr


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   call_dmr/index.rst
   call_dms/index.rst
   parse_methylpy/index.rst
   rms_test/index.rst


Package Contents
----------------

.. py:function:: call_dmr(output_dir, replicate_label=None, p_value_cutoff=0.001, frac_delta_cutoff=0.2, max_dist=250, residual_quantile=0.6, corr_cutoff=0.3, dms_ratio=0.8, cpu=1, chrom=None)

   Call DMR from DMS results.

   :param output_dir:
   :param replicate_label:
   :param p_value_cutoff:
   :param frac_delta_cutoff:
   :param max_dist:
   :param residual_quantile:
   :param corr_cutoff:
   :param dms_ratio:
   :param cpu:
   :param chrom:


.. py:function:: collapse_replicates(region_ds, replicate_label, state_da='dmr_state')


.. py:function:: call_dms(output_dir, allc_paths, samples, chrom_size_path, cpu=1, max_row_count=50, n_permute=3000, min_pvalue=0.01, region=None)

   Call DMS from multiple ALLC files

   :param output_dir:
   :param allc_paths:
   :param samples:
   :param chrom_size_path:
   :param cpu:
   :param max_row_count:
   :param n_permute:
   :param min_pvalue:
   :param region:


.. py:function:: methylpy_to_region_ds(dmr_path, output_dir)


