:py:mod:`ALLCools._doc`
=======================

.. py:module:: ALLCools._doc


Module Contents
---------------

.. py:data:: idx_doc
   :annotation: = If true, save an methylpy chromosome index for back compatibility. If you only use methylpy to...

   

.. py:data:: allc_path_doc
   :annotation: = Path to 1 ALLC file

   

.. py:data:: allc_paths_doc
   :annotation: = Single ALLC path contain wildcard OR multiple space separated ALLC paths OR a file contains 1...

   

.. py:data:: allc_table_doc
   :annotation: = Contain all the ALLC file information in two tab-separated columns: 1. file_uid, 2. file_path. No header

   

.. py:data:: binarize_doc
   :annotation: = If set, binarize each single site in each individual ALLC file. This means each cytosine will...

   

.. py:data:: bin_sizes_doc
   :annotation: = Fix-size genomic bins can be defined by bin_sizes and chrom_size_path. Space separated sizes of...

   

.. py:data:: bw_bin_sizes_doc
   :annotation: = Bin size of the BigWig files.

   

.. py:data:: chrom_size_path_doc
   :annotation: = Path to UCSC chrom size file. This can be generated from the genome fasta or downloaded via UCSC...

   

.. py:data:: compress_level_doc
   :annotation: = Compression level for the output file

   

.. py:data:: cov_cutoff_doc
   :annotation: = Max cov filter for a single site in ALLC. Sites with cov > cov_cutoff will be skipped.

   

.. py:data:: cpu_basic_doc
   :annotation: = Number of processes to use in parallel.

   

.. py:data:: mc_contexts_doc
   :annotation: = Space separated mC context patterns to extract from ALLC. The context length should be the same...

   

.. py:data:: mc_context_mcad_doc
   :annotation: = mC context pattern to extract from ALLC. Context pattern follows IUPAC nucleotide code, e.g. N...

   

.. py:data:: reference_fasta_doc
   :annotation: = Path to 1 genome reference FASTA file (the one used for mapping), use samtools fadix to build...

   

.. py:data:: region_bed_names_doc
   :annotation: = Space separated names for each BED file provided in region_bed_paths.

   

.. py:data:: region_bed_paths_doc
   :annotation: = Arbitrary genomic regions can be defined in several BED files to count on. Space separated paths...

   

.. py:data:: region_bed_path_mcad_doc
   :annotation: = Arbitrary genomic regions can be defined in one BED file to count on. The fourth column of the...

   

.. py:data:: region_doc
   :annotation: = Only extract records from certain genome region(s) via tabix, multiple region can be provided in...

   

.. py:data:: remove_additional_chrom_doc
   :annotation: = Whether to remove rows with unknown chromosome instead of raising KeyError

   

.. py:data:: rna_table_doc
   :annotation: = This is only for mCT data when we have RNA BAM file for each single cell. Contain all the RNA...

   

.. py:data:: snp_doc
   :annotation: = If true, means the input allc contain snp information, and the allc processing will take care that.

   

.. py:data:: split_strand_doc
   :annotation: = If true, Watson (+) and Crick (-) strands will be count separately

   

.. py:data:: strandness_doc
   :annotation: = What to do with strand information, possible values are: 1. both: save +/- strand together in...

   

.. py:data:: generate_dataset_doc
   :annotation: = Generate MCDS dataset from a list of ALLC files (recorded in the allc_table). Multiple region...

   

.. py:data:: generate_dataset_obs_dim_doc
   :annotation: = Name of the observation dimension.

   

.. py:data:: generate_dataset_chunk_size_doc
   :annotation: = Chunk allc_table with chunk_size when generate dataset in parallel

   

.. py:data:: generate_dataset_regions_doc
   :annotation: = Definition of genomic regions in the form of "--regions {region_name} {region_definition}". This...

   

.. py:data:: generate_dataset_quantifiers_doc
   :annotation: = Definition of genome region quantifiers in the form of "--quantifiers {region_name} {quant_type}...

   

.. py:data:: table_to_allc_doc
   :annotation: = Convert different kinds of methylation table into ALLC format. Currently, only plain text table...

   

.. py:data:: table_to_allc_input_path
   :annotation: = input path of the table

   

.. py:data:: table_to_allc_output_prefix
   :annotation: = output prefix of the ALLC table

   

.. py:data:: table_to_allc_sep
   :annotation: = character to separate columns in the table

   

.. py:data:: table_to_allc_header
   :annotation: = Whether the table contains header line or not

   

.. py:data:: table_to_allc_chunk_size
   :annotation: = chunk_size to perform conversion

   

.. py:data:: table_to_allc_chrom
   :annotation: = the chromosome column number, 0-based index

   

.. py:data:: table_to_allc_pos
   :annotation: = the position column number, 0-based index

   

.. py:data:: table_to_allc_strand
   :annotation: = the strand column number, 0-based index. If not provided, will infer automatically based on the...

   

.. py:data:: table_to_allc_context
   :annotation: = the cytosine context column number, 0-based index. If not provided, will inter automatically...

   

.. py:data:: table_to_allc_mc
   :annotation: = the methylated cytosine count column number, 0-based index.

   

.. py:data:: table_to_allc_uc
   :annotation: = the unmethylated cytosine count column number, 0-based index.

   

.. py:data:: table_to_allc_cov
   :annotation: = the total cytosine coverage count column number, 0-based index.

   

.. py:data:: table_to_allc_mc_frac
   :annotation: = the methylation fraction column number, 0-based index.

   

.. py:data:: table_to_allc_pseudo_count
   :annotation: = Use this pseudo_count number as the total cytosine coverage count, if the "cov" column is...

   

.. py:data:: table_to_allc_fasta_path
   :annotation: = the genome FASTA file path, required if either "strand" or "context" column is missing.

   

.. py:data:: table_to_allc_num_upstream_bases
   :annotation: = number of up stream bases to include when get cytosine context.

   

.. py:data:: table_to_allc_num_downstream_bases
   :annotation: = number of down stream bases to include when get cytosine context.

   

.. py:data:: table_to_allc_add_chr
   :annotation: = whether add "chr" before the chromosome name.

   

.. py:data:: table_to_allc_sort
   :annotation: = whether sort the ALLC table after conversion.

   

.. py:function:: doc_params(**kwds)

   Docstrings should start with "" in the first line for proper formatting.



