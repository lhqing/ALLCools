:py:mod:`ALLCools.dmr.call_dmr`
===============================

.. py:module:: ALLCools.dmr.call_dmr


Module Contents
---------------

.. py:function:: _call_dmr_single_chrom(dms_dir, output_dir, chrom, p_value_cutoff=0.001, frac_delta_cutoff=0.2, max_dist=250, residual_quantile=0.6, corr_cutoff=0.3, dms_ratio=0.8)

   Call DMR for single chromosome, see call_dmr for doc


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


