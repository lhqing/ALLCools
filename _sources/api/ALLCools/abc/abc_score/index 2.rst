:py:mod:`ALLCools.abc.abc_score`
================================

.. py:module:: ALLCools.abc.abc_score


Module Contents
---------------

.. py:function:: activity_score_atac_mcg(all_values)


.. py:function:: scan_single_bw(bw_path, bed, value_type='sum')


.. py:function:: preprocess_bed(input_bed, blacklist_path, chrom_sizes, temp_dir)


.. py:class:: ABCModel(cool_url, enhancer_bed_path, tss_bed_path, output_prefix, epi_mark: dict, calculation_mode, blacklist_path=None, enhancer_size=500, promoter_size=500, max_dist=5000000, min_score_cutoff=0.02, balance=False, cleanup=True, cpu=1)

   .. py:method:: _enhancer_bed(self)


   .. py:method:: _tss_bed(self)


   .. py:method:: _tss_bin(self)


   .. py:method:: _promoter_bed(self)


   .. py:method:: _promoter_and_enhancer_bed(self)


   .. py:method:: _standard_bed(self, bed, name, slop)


   .. py:method:: _single_chrom_activity_score(self, chrom_pe_bed)


   .. py:method:: _single_chrom_abc_score(self, chrom_mat, chrom_tss_bin, chrom_pe_bed)


   .. py:method:: calculate(self, balance=False, cpu=1)



