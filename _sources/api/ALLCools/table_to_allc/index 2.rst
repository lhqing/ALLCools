:py:mod:`ALLCools.table_to_allc`
================================

.. py:module:: ALLCools.table_to_allc


Module Contents
---------------

.. py:function:: mode_mc_cov(table, mc, cov)


.. py:function:: mode_mc_uc(table, mc, uc)


.. py:function:: mode_uc_cov(table, uc, cov)


.. py:function:: mode_mc_frac_cov(table, mc_frac, cov)


.. py:function:: mode_mc_frac_mc(table, mc_frac, mc)


.. py:function:: mode_mc_frac_uc(table, mc_frac, uc)


.. py:function:: mode_mc_frac_pseudo_count(table, mc_frac, pseudo_count)


.. py:function:: get_strand(chrom, pos, fasta)


.. py:function:: get_context(chrom, pos, strand, fasta, num_upstream_bases, num_downstream_bases)


.. py:function:: get_strand_and_context(table, chrom, pos, strand, context, fasta_path, num_upstream_bases, num_downstream_bases)


.. py:function:: dataframe_to_allc(table, add_chr=False, chrom=0, pos=1, strand=None, context=None, fasta_path=None, chrom_sizes=None, mc=None, uc=None, cov=None, mc_frac=None, pseudo_count=1, num_upstream_bases=0, num_downstream_bases=2)


.. py:function:: table_to_allc(input_path, output_prefix, sep='\t', header=None, chunk_size=100000, chrom=0, pos=1, strand=None, context=None, mc=None, uc=None, cov=None, mc_frac=None, pseudo_count=1, fasta_path=None, num_upstream_bases=0, num_downstream_bases=2, add_chr=False, sort=True)

   {table_to_allc_doc}

   :param input_path: {table_to_allc_input_path}
   :param output_prefix: {table_to_allc_output_prefix}
   :param sep: {table_to_allc_sep}
   :param header: {table_to_allc_header}
   :param chunk_size: {table_to_allc_chunk_size}
   :param chrom: {table_to_allc_chrom}
   :param pos: {table_to_allc_pos}
   :param strand: {table_to_allc_strand}
   :param context: {table_to_allc_context}
   :param mc: {table_to_allc_mc}
   :param uc: {table_to_allc_uc}
   :param cov: {table_to_allc_cov}
   :param mc_frac: {table_to_allc_mc_frac}
   :param pseudo_count: {table_to_allc_pseudo_count}
   :param fasta_path: {table_to_allc_fasta_path}
   :param num_upstream_bases: {table_to_allc_num_upstream_bases}
   :param num_downstream_bases: {table_to_allc_num_downstream_bases}
   :param add_chr: {table_to_allc_add_chr}
   :param sort: {table_to_allc_sort}

   :returns:
   :rtype: output path of the converted ALLC file.


