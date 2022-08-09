:py:mod:`ALLCools._merge_allc`
==============================

.. py:module:: ALLCools._merge_allc

.. autoapi-nested-parse::

   Some of the functions are modified from methylpy https://github.com/yupenghe/methylpy

   Original Author: Yupeng He



Module Contents
---------------

.. py:data:: log
   

   

.. py:data:: DEFAULT_MAX_ALLC
   :annotation: = 150

   

.. py:data:: PROCESS
   

   

.. py:class:: _ALLC(path, region)

   .. py:method:: readline(self)


   .. py:method:: close(self)



.. py:function:: _increase_soft_fd_limit()

   Increase soft file descriptor limit to hard limit,
   this is the maximum a process can do
   Use this in merge_allc, because for single cell, a lot of file need to be opened

   Some useful discussion
   https://unix.stackexchange.com/questions/36841/why-is-number-of-open-files-limited-in-linux
   https://docs.python.org/3.6/library/resource.html
   https://stackoverflow.com/questions/6774724/why-python-has-limit-for-count-of-file-handles/6776345


.. py:function:: _batch_merge_allc_files_tabix(allc_files, out_file, chrom_size_file, bin_length, cpu=10, binarize=False, snp=False)


.. py:function:: _merge_allc_files_tabix(allc_files, out_file, chrom_size_file, query_region=None, buffer_line_number=10000, binarize=False)


.. py:function:: _merge_allc_files_tabix_with_snp_info(allc_files, out_file, chrom_size_file, query_region=None, buffer_line_number=10000, binarize=False)


.. py:function:: merge_allc_files(allc_paths, output_path, chrom_size_path, bin_length=10000000, cpu=10, binarize=False, snp=False)

   Merge N ALLC files into 1 ALLC file.

   :param allc_paths: {allc_paths_doc}
   :param output_path: Path to the output merged ALLC file.
   :param chrom_size_path: {chrom_size_path_doc}
   :param bin_length: Length of the genome bin in each parallel job, large number means more memory usage.
   :param cpu: {cpu_basic_doc}
               The real CPU usage is ~1.5 times than this number,
               due to the sub processes of handling ALLC files using tabix/bgzip.
               Monitor the CPU and Memory usage when running this function.
   :param binarize: {binarize_doc}
   :param snp: {snp_doc}


