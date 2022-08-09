:py:mod:`ALLCools.gtf.Gtf`
==========================

.. py:module:: ALLCools.gtf.Gtf


Module Contents
---------------

.. py:function:: create_gtf_db(gtf_path, disable_infer_genes=True, disable_infer_transcripts=True)


.. py:class:: Gtf(*args, **kwargs)

   Bases: :py:obj:`gffutils.FeatureDB`

   .. py:method:: get_gene_name(self, feature_id)


   .. py:method:: _select_longest_id(self, id_list)


   .. py:method:: get_gene_id_by_name(self, gene_name, select_longest=True)


   .. py:method:: _convert_to_id(self, name)


   .. py:method:: get_gene_feature(self, gene)


   .. py:method:: get_gene_length(self, gene)


   .. py:method:: get_gene_transcripts(self, gene, featuretype='transcript')


   .. py:method:: get_gene_exons(self, gene, featuretype='exon')


   .. py:method:: get_gene_tss(self, gene, transcript_featuretype='transcript')


   .. py:method:: get_gene_promoter(self, gene, chrom_sizes_path, slop=500, transcript_featuretype='transcript')


   .. py:method:: subset_db(self, output_path, genes=None, region_bed=None, span=500000, disable_infer_genes=True, disable_infer_transcripts=True)



