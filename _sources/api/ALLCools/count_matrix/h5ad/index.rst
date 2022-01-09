:py:mod:`ALLCools.count_matrix.h5ad`
====================================

.. py:module:: ALLCools.count_matrix.h5ad


Module Contents
---------------

.. py:function:: _bin_count_table_to_csr_npz(bin_count_tables, bin_size, chrom_size_path, output_prefix, compression=True)


.. py:function:: _csr_matrix_to_anndata(matrix_paths, output_path, obs_names, chrom_size_path, bin_size, mc_type, count_type, step_size, strandness, compression=None, compression_opts=None)


.. py:function:: aggregate_region_count_to_paired_anndata(count_tables, output_prefix, chrom_size_path, bin_size, mc_type, count_type, strandness, compression='gzip', file_uids=None, max_obj=3072, cpu=3)


.. py:function:: _transform_single_h5ad(adata_path, output_path, chrom_size_path, bin_size, step_size, window_size, compression)

   Resize non-overlap chrom bin count adata


