:py:mod:`ALLCools.sandbox.motif.fimo`
=====================================

.. py:module:: ALLCools.sandbox.motif.fimo

.. autoapi-nested-parse::

   Run FIMO to scan motif over FASTA sequences



Module Contents
---------------

.. py:function:: _fimo_runner(motif_path, fasta_path, output_path, path_to_fimo='', raw_score_thresh=6, raw_p_value_thresh=0.0005, is_genome_fasta=False, top_n=300000)

   Run fimo for a single motif over single fasta file


.. py:function:: _scan_motif_over_fasta(meme_motif_file, fasta_path, cpu, output_dir, path_to_fimo='', raw_score_thresh=7, raw_p_value_thresh=0.0005, top_n=300000, is_genome_fasta=False)


.. py:function:: _aggregate_motif_beds(bed_file_paths, output_path, cpu=1, sort_mem_gbs=1)


