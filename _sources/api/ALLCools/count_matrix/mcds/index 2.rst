:py:mod:`ALLCools.count_matrix.mcds`
====================================

.. py:module:: ALLCools.count_matrix.mcds


Module Contents
---------------

.. py:data:: DEFAULT_MCDS_DTYPE
   

   

.. py:function:: clip_too_large_cov(data_df, machine_max)


.. py:function:: _region_count_table_to_csr_npz(region_count_tables, region_id_map, output_prefix, compression=True, dtype=DEFAULT_MCDS_DTYPE)

   helper func of _aggregate_region_count_to_mcds

   Take a list of region count table paths, read,
   aggregate them into a 2D sparse matrix and save the mC and COV separately.
   This function don't take care of any path selection, but assume all region_count_table is homogeneous type
   It return the saved file path


.. py:function:: _csr_matrix_to_dataarray(matrix_table, row_name, row_index, col_name, col_index, other_dim_info)

   helper func of _aggregate_region_count_to_mcds

   This function aggregate sparse array files into a single xarray.DataArray,
   combining cell chunks, mc/cov count type together.
   The matrix_table provide all file paths, each row is for a cell chunk, with mc and cov matrix path separately.


.. py:function:: _aggregate_region_count_to_mcds(output_dir, dataset_name, chunk_size=100, row_name='cell', cpu=1, dtype=DEFAULT_MCDS_DTYPE)

   This function aggregate all the region count table into a single mcds


.. py:function:: generate_mcds(allc_table, output_prefix, chrom_size_path, mc_contexts, rna_table=None, split_strand=False, bin_sizes=None, region_bed_paths=None, region_bed_names=None, cov_cutoff=9999, cpu=1, remove_tmp=True, max_per_mcds=3072, cell_chunk_size=100, dtype=DEFAULT_MCDS_DTYPE, binarize=False, engine='zarr')

   Generate MCDS from a list of ALLC file provided with file id.

   :param allc_table: {allc_table_doc}
   :param output_prefix: Output prefix of the MCDS
   :param chrom_size_path: {chrom_size_path_doc}
   :param mc_contexts: {mc_contexts_doc}
   :param rna_table: {rna_table_doc}
   :param split_strand: {split_strand_doc}
   :param bin_sizes: {bin_sizes_doc}
   :param region_bed_paths: {region_bed_paths_doc}
   :param region_bed_names: {region_bed_names_doc}
   :param cov_cutoff: {cov_cutoff_doc}
   :param cpu: {cpu_basic_doc}
   :param remove_tmp: Whether to remove the temp directory for generating MCDS
   :param max_per_mcds: Maximum number of ALLC files to aggregate into 1 MCDS, if number of ALLC provided > max_per_mcds,
                        will generate MCDS in chunks, with same prefix provided.
   :param cell_chunk_size: Size of cell chunk in parallel aggregation. Do not have any effect on results.
                           Large chunksize needs large memory.
   :param dtype: Data type of MCDS count matrix. Default is np.uint32.
                 For single cell feature count, this can be set to np.uint16, which means the value is 0-65536.
                 The values exceed max will be clipped.
   :param binarize: {binarize_doc}
   :param engine: use zarr or netcdf to store dataset, default is zarr


