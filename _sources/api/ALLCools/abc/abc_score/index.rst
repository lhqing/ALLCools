:py:mod:`ALLCools.abc.abc_score`
================================

.. py:module:: ALLCools.abc.abc_score


Module Contents
---------------

.. py:function:: _activity_score_atac_mcg(all_values)


.. py:function:: _scan_single_bw(bw_path, bed, value_type='sum')


.. py:function:: _preprocess_bed(input_bed, blacklist_path, chrom_sizes, temp_dir)


.. py:class:: ABCModel(cool_url, enhancer_bed_path, tss_bed_path, output_prefix, epi_mark: dict, calculation_mode, blacklist_path=None, enhancer_size=500, promoter_size=500, max_dist=5000000, min_score_cutoff=0.02, balance=False, cleanup=True, cpu=1)

   Calculate ABC score between promoter and enhancer.

   .. py:method:: _enhancer_bed()


   .. py:method:: _tss_bed()


   .. py:method:: _tss_bin()


   .. py:method:: _promoter_bed()


   .. py:method:: _promoter_and_enhancer_bed()


   .. py:method:: _standard_bed(bed, name, slop)


   .. py:method:: _single_chrom_activity_score(chrom_pe_bed)


   .. py:method:: _single_chrom_abc_score(chrom_mat, chrom_tss_bin, chrom_pe_bed)


   .. py:method:: calculate(balance=False, cpu=1)

      Calculate the ABC score.



