:py:mod:`ALLCools._bam_to_allc`
===============================

.. py:module:: ALLCools._bam_to_allc

.. autoapi-nested-parse::

   This file is modified from methylpy https://github.com/yupenghe/methylpy.

   Author: Yupeng He

   Notes added on 07/03/2022 - Difference between bismark and hisat-3n
   Read mpileup doc first: http://www.htslib.org/doc/samtools-mpileup.html

   For Bismark SE mapping:
   Bismark converted the reads orientation based on their conversion type.
   C to T conversion reads are always map to the forward strand, regardless of the R1 R2
   G to A conversion reads are always map to the reverse strand, regardless of the R1 R2
   Therefore, in the original bam-to-allc function, we simply consider the strandness in
   the pileup format to distinguish which C needs to be counted.
   There are two situations:
   1. If the ref base is C, the read need to map to forward strand in order to be counted
   [.ATCGN] corresponding to the forward strand
   2. If the ref base is G, the read need to map to reverse strand in order to be counted
   [,atcgn] corresponding to the reverse strand

   For Hisat-3n PE or Biskarp PE mapping:
   R1 R2 are mapped as their original orientation, therefore,
   both the C to T and G to A conversion reads can have forward and reverse strand alignment.
   We can not distinguish the conversion type by strand in input bam.
   Here I add a check using the YZ tag of hisat-3n BAM file or XG tag of bismark BAM file.
   If the YZ tag is "+" or XG tag is "CT", the read is C to T conversion, I change the flag to forward mapping
   no matter R1 or R2 by read.is_forward = True
   If the YZ tag is "-" or XG tag is "GA", the read is G to A conversion, I change the flag to reverse mapping
   no matter R1 or R2 by read.is_forward = False
   In this case, the read orientation is the same as bismark bam file, and the following base
   count code no need to change.



Module Contents
---------------

.. py:data:: log
   

   

.. py:function:: _is_read_ct_conversion_hisat3n(read)


.. py:function:: _is_read_ct_conversion_bismark(read)


.. py:function:: _convert_bam_strandness(in_bam_path, out_bam_path)


.. py:function:: _read_faidx(faidx_path)

   Read fadix of reference fasta file.

   samtools fadix ref.fa


.. py:function:: _get_chromosome_sequence_upper(fasta_path, fai_df, query_chrom)

   Read a whole chromosome sequence into memory.


.. py:function:: _get_bam_chrom_index(bam_path)


.. py:function:: _bam_to_allc_worker(bam_path, reference_fasta, fai_df, output_path, region=None, num_upstr_bases=0, num_downstr_bases=2, buffer_line_number=100000, min_mapq=0, min_base_quality=1, compress_level=5, tabix=True, save_count_df=False)

   None parallel bam_to_allc worker function, call by bam_to_allc.


.. py:function:: _aggregate_count_df(count_dfs)


.. py:function:: bam_to_allc(bam_path, reference_fasta, output_path=None, cpu=1, num_upstr_bases=0, num_downstr_bases=2, min_mapq=10, min_base_quality=20, compress_level=5, save_count_df=False, convert_bam_strandness=False)

   Generate 1 ALLC file from 1 position sorted BAM file via samtools mpileup.

   :param bam_path: Path to 1 position sorted BAM file
   :param reference_fasta: {reference_fasta_doc}
   :param output_path: Path to 1 output ALLC file
   :param cpu: {cpu_basic_doc} DO NOT use cpu > 1 for single cell ALLC generation.
               Parallel on cell level is better for single cell project.
   :param num_upstr_bases: Number of upstream base(s) of the C base to include in ALLC context column,
                           usually use 0 for normal BS-seq, 1 for NOMe-seq.
   :param num_downstr_bases: Number of downstream base(s) of the C base to include in ALLC context column,
                             usually use 2 for both BS-seq and NOMe-seq.
   :param min_mapq: Minimum MAPQ for a read being considered, samtools mpileup parameter, see samtools documentation.
   :param min_base_quality: Minimum base quality for a base being considered, samtools mpileup parameter,
                            see samtools documentation.
   :param compress_level: {compress_level_doc}
   :param save_count_df: If true, save an ALLC context count table next to ALLC file.
   :param convert_bam_strandness: {convert_bam_strandness_doc}

   :returns: a pandas.DataFrame for overall mC and cov count separated by mC context.
   :rtype: count_df


