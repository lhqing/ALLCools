:py:mod:`ALLCools.motif`
========================

.. py:module:: ALLCools.motif


Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   default_motif_set/index.rst


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   motifs/index.rst
   parse_meme/index.rst
   snakemake/index.rst


Package Contents
----------------

.. py:class:: MotifSet(motifs, meta_table=None, motif_cluster_col='cluster_id')

   .. py:method:: calculate_threshold(method='balanced', cpu=1, threshold_value=1000)


   .. py:method:: _multi_seq_motif_scores(seq_list, save_path, region_dim, motif_dim, dtype)


   .. py:method:: _run_motif_scan_chunks(fasta_path, output_dir, region_dim, motif_dim, dtype, chunk_size=10000, cpu=1)


   .. py:method:: _aggregate_motif_clusters(output_dir, motif_dim, region_dim, dtype, chunk_size=1000000)


   .. py:method:: scan_motifs(fasta_path, output_dir, cpu, region_dim, motif_dim='motif', combine_cluster=True, dtype='uint16', chunk_size=10000)



.. py:function:: get_default_motif_set(database='three_databases')

   Get motif set from pre-defined motif database.


