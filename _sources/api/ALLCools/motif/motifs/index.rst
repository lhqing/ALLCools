:py:mod:`ALLCools.motif.motifs`
===============================

.. py:module:: ALLCools.motif.motifs


Module Contents
---------------

.. py:function:: _get_motif_threshold(motif, method, threshold_value)


.. py:function:: _single_seq_single_motif_max_score(seq, motif, threshold)


.. py:class:: MotifSet(motifs, meta_table=None, motif_cluster_col='cluster_id')

   .. py:method:: calculate_threshold(method='balanced', cpu=1, threshold_value=1000)


   .. py:method:: _multi_seq_motif_scores(seq_list, save_path, region_dim, motif_dim, dtype)


   .. py:method:: _run_motif_scan_chunks(fasta_path, output_dir, region_dim, motif_dim, dtype, chunk_size=10000, cpu=1)


   .. py:method:: _aggregate_motif_clusters(output_dir, motif_dim, region_dim, dtype, chunk_size=1000000)


   .. py:method:: scan_motifs(fasta_path, output_dir, cpu, region_dim, motif_dim='motif', combine_cluster=True, dtype='uint16', chunk_size=10000)



.. py:class:: Motif(alphabet='ACGT', instances=None, counts=None)

   Bases: :py:obj:`Bio.motifs.Motif`

   A class representing sequence motifs.


