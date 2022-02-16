:py:mod:`ALLCools.sandbox.motif.ame`
====================================

.. py:module:: ALLCools.sandbox.motif.ame


Module Contents
---------------

.. py:function:: _single_ame_runner(fasta_path, motif_paths, ame_kws_str='')


.. py:function:: ame(bed_file, motif_files, output_path, reference_fasta, background_files=None, cpu=1, standard_length=None, chrom_size_path=None, ame_kws=None, sample=None, seed=1)

   Motif enrichment analysis with AME from MEME Suite.

   :param bed_file: Single bed file input to calculate motif enrichment
   :param motif_files: One or multiple motif files in MEME format
   :param output_path: Output path of the AME result
   :param reference_fasta: Reference fasta of the bed_file and background_files
   :param background_files: One or multiple bed files contain region as the control set of motif enrichment
   :param cpu: Number of CPUs to parallel AME
   :param standard_length: If not None, will standard the input region and control region (if provided) in to same standard_length,
                           each standardized region is centered to the original region center.
   :param chrom_size_path: {chrom_size_path_doc}
                           Must provide is standard_length is not None.
   :param ame_kws: Additional options that will pass to AME, provide either a dict or string.
                   If string provided, will directly insert into the final AME command,
                   see AME documentation about AME options.
                   e.g. '--seed 1 --method fisher --scoring avg'
   :param sample: Sample input regions to this number to reduce run time or test
   :param seed: Seed for random sample input regions, only apply when sample is not None


