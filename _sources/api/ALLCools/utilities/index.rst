:py:mod:`ALLCools.utilities`
============================

.. py:module:: ALLCools.utilities


Module Contents
---------------

.. py:data:: IUPAC_TABLE
   

   

.. py:data:: COMPLIMENT_BASE
   

   

.. py:function:: reverse_compliment(seq)


.. py:function:: get_allc_chroms(allc_path)


.. py:function:: parse_mc_pattern(pattern: str) -> set

   parse mC context pattern


.. py:function:: parse_chrom_size(path, remove_chr_list=None)

   support simple UCSC chrom size file, or .fai format (1st and 2nd columns same as chrom size file)

   return chrom:length dict


.. py:function:: chrom_dict_to_id_index(chrom_dict, bin_size)


.. py:function:: get_bin_id(chrom, chrom_index_dict, bin_start, bin_size) -> int


.. py:function:: genome_region_chunks(chrom_size_path: str, bin_length: int = 10000000, combine_small: bool = True) -> List[str]

   Split the whole genome into bins, where each bin is {bin_length} bp. Used for tabix region query

   :param chrom_size_path: Path of UCSC genome size file
   :param bin_length: length of each bin
   :param combine_small: whether combine small regions into one record

   :returns:
   :rtype: list of records in tabix query format


.. py:function:: parse_file_paths(input_file_paths: Union[str, list]) -> list


.. py:function:: get_md5(file_path)


.. py:function:: check_tbi_chroms(file_path, genome_dict, same_order=False)


.. py:function:: generate_chrom_bin_bed_dataframe(chrom_size_path: str, window_size: int, step_size: int = None) -> pandas.DataFrame

   Generate BED format dataframe based on UCSC chrom size file and window_size
   return dataframe contain 3 columns: chrom, start, end. The index is 0 based continue bin index.


.. py:function:: profile_allc(allc_path, drop_n=True, n_rows=1000000, output_path=None)

   Generate some summary statistics of 1 ALLC.
   1e8 rows finish in about 5 min.

   :param allc_path: {allc_path_doc}
   :param drop_n: Whether to drop context that contain N, such as CCN.
                  This is usually very rare and need to be dropped.
   :param n_rows: Number of rows to calculate the profile from.
                  The default number is usually sufficient to get pretty precise assumption.
   :param output_path: Path of the output file. If None, will save the profile next to input ALLC file.


.. py:function:: is_gz_file(filepath)

   Check if a file is gzip file, bgzip also return True
   Learnt from here: https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed


.. py:function:: tabix_allc(allc_path, reindex=False)

   a simple wrapper of tabix command to index 1 ALLC file

   :param allc_path: {allc_path_doc}
   :param reindex: If True, will force regenerate the ALLC index.


.. py:function:: standardize_allc(allc_path, chrom_size_path, compress_level=5, remove_additional_chrom=False)

   Standardize 1 ALLC file by checking:
       1. No header in the ALLC file;
       2. Chromosome names in ALLC must be same as those in the chrom_size_path file, including "chr";
       3. Output file will be bgzipped with .tbi index
       4. Remove additional chromosome (remove_additional_chrom=True) or
          raise KeyError if unknown chromosome found (default)

   :param allc_path: {allc_path_doc}
   :param chrom_size_path: {chrom_size_path_doc}
   :param compress_level: {compress_level_doc}
   :param remove_additional_chrom: {remove_additional_chrom_doc}


.. py:function:: _transfer_bin_size(bin_size: int) -> str

   Get proper str for a large bin_size


.. py:function:: parse_dtype(dtype)


.. py:function:: binary_count(mc, cov)


