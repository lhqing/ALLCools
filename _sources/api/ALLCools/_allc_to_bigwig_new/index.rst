:py:mod:`ALLCools._allc_to_bigwig_new`
======================================

.. py:module:: ALLCools._allc_to_bigwig_new


Module Contents
---------------

.. py:class:: ContextCounter(mc_contexts)

   .. py:method:: add(self, context, mc, cov)



.. py:class:: StrandContextCounter(mc_contexts)

   .. py:method:: add(self, context, strand, mc, cov)



.. py:function:: write_entry(counter, context_handle, mc_contexts, strandness, chrom, bin_start, bin_size)


.. py:function:: allc_to_bigwig(allc_path, output_prefix, bin_size, mc_contexts, chrom_size_path, strandness)

   Generate BigWig files from one ALLC file.

   :param allc_path: {allc_path_doc}
   :param output_prefix: Path prefix of the output BigWig file.
   :param bin_size: {bw_bin_sizes_doc}
   :param mc_contexts: {mc_contexts_doc}
   :param strandness: {strandness_doc}
   :param chrom_size_path: {chrom_size_path_doc}
                           If chrom_size_path provided, will use it to extract ALLC with chrom order,
                           but if region provided, will ignore this.


