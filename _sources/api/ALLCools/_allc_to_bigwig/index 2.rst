:py:mod:`ALLCools._allc_to_bigwig`
==================================

.. py:module:: ALLCools._allc_to_bigwig

.. autoapi-nested-parse::

   This file is modified from methylpy https://github.com/yupenghe/methylpy



Module Contents
---------------

.. py:data:: log
   

   

.. py:function:: _allc_to_bedgraph(allc_path, out_prefix, chrom_size_path, remove_additional_chrom=False, bin_size=50)

   Simply calculate cov and mc_rate for fixed genome bins. No mC context filter.


.. py:function:: _bedgraph_to_bigwig(input_file, chrom_size_path, path_to_wigtobigwig, remove_bedgraph=True)


.. py:function:: allc_to_bigwig(allc_path, output_prefix, chrom_size_path, mc_contexts, split_strand=False, bin_size=50, remove_additional_chrom=False, region=None, cov_cutoff=9999, path_to_wigtobigwig='', remove_temp=True, cpu=1)

   Generate bigwig file(s) from 1 ALLC file.

   :param allc_path: {allc_path_doc}
   :param output_prefix: Output prefix of the bigwig file(s)
   :param chrom_size_path: {chrom_size_path_doc}
   :param mc_contexts: {mc_contexts_doc}
   :param split_strand: {split_strand_doc}
   :param bin_size: Minimum bin size of bigwig file
   :param remove_additional_chrom: {remove_additional_chrom_doc}
   :param region: {region_doc}
   :param cov_cutoff: {cov_cutoff_doc}
   :param path_to_wigtobigwig: Path to wigtobigwig to allow allcools to find it
   :param remove_temp: debug parameter, whether to remove the temp file or not
   :param cpu: Number of cores to use


