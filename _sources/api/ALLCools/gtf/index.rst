:py:mod:`ALLCools.gtf`
======================

.. py:module:: ALLCools.gtf


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   Gtf/index.rst
   utilities/index.rst


Package Contents
----------------

.. py:class:: Gtf(*args, **kwargs)

   Bases: :py:obj:`gffutils.FeatureDB`

   .. py:method:: get_gene_name(feature_id)


   .. py:method:: _select_longest_id(id_list)


   .. py:method:: get_gene_id_by_name(gene_name, select_longest=True)


   .. py:method:: _convert_to_id(name)


   .. py:method:: get_gene_feature(gene)


   .. py:method:: get_gene_length(gene)


   .. py:method:: get_gene_transcripts(gene, featuretype='transcript')


   .. py:method:: get_gene_exons(gene, featuretype='exon')


   .. py:method:: get_gene_tss(gene, transcript_featuretype='transcript')


   .. py:method:: get_gene_promoter(gene, chrom_sizes_path, slop=500, transcript_featuretype='transcript')


   .. py:method:: subset_db(output_path, genes=None, region_bed=None, span=500000, disable_infer_genes=True, disable_infer_transcripts=True)



.. py:function:: create_gtf_db(gtf_path, disable_infer_genes=True, disable_infer_transcripts=True)


.. py:function:: read_gtf(gtf_path)

   Read GTF file.


.. py:function:: subset_gtf(gtf, regions, output_path=None, select_feature=None)

   Subset GTF file by genomic regions.


