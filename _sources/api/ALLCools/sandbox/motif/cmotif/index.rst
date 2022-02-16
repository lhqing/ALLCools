:py:mod:`ALLCools.sandbox.motif.cmotif`
=======================================

.. py:module:: ALLCools.sandbox.motif.cmotif


Module Contents
---------------

.. py:function:: _generate_c_motif_database(motif_bed, fasta_path, output_dir, bin_size=10000000)


.. py:function:: generate_cmotif_database(bed_file_paths, reference_fasta, motif_files, output_dir, slop_b=None, chrom_size_path=None, cpu=1, sort_mem_gbs=5, path_to_fimo='', raw_score_thresh=8.0, raw_p_value_thresh=0.0002, top_n=300000, cmotif_bin_size=10000000)

   Generate lookup table for motifs all the cytosines belongs to.
   BED files are used to limit cytosine scan in certain regions.
   Scanning motif over whole genome is very noisy, better scan it in some functional part of genome.
   The result files will be in the output

   :param bed_file_paths: Paths of bed files. Multiple bed will be merged to get a final region set.
                          The motif scan will only happen on the regions defined in these bed files.
   :param reference_fasta: FASTA file path of the genome to scan
   :param motif_files: MEME motif files that contains all the motif information.
   :param output_dir: Output directory of C-Motif database
   :param slop_b: Whether add slop to both ends of bed files.
   :param chrom_size_path: {chrom_size_path_doc}
                           Needed if slop_b is not None
   :param cpu: {cpu_basic_doc}
   :param sort_mem_gbs: Maximum memory usage in GBs when sort bed files
   :param path_to_fimo: Path to fimo executable, if fimo is not in PATH
   :param raw_score_thresh: Threshold of raw motif match likelihood score, see fimo doc for more info.
   :param raw_p_value_thresh: Threshold of raw motif match P-value, see fimo doc for more info.
   :param top_n: If too much motif found, will order them by likelihood score and keep top matches.
   :param cmotif_bin_size: Bin size of single file in C-Motif database. No impact on results, better keep the default.


