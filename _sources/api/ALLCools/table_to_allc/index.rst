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


