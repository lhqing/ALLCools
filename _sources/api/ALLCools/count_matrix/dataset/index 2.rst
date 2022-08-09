:py:mod:`ALLCools.count_matrix.dataset`
=======================================

.. py:module:: ALLCools.count_matrix.dataset


Module Contents
---------------

.. py:data:: ALLOW_QUANT_TYPES
   :annotation: = ['count', 'hypo-score', 'hyper-score']

   

.. py:function:: bin_sf(cov, mc, p)


.. py:function:: cell_sf(cell_count_df)


.. py:class:: Quant(mc_types, quant_type, kwargs)


.. py:class:: CountQuantifier(mc_types)

   Sum count by mC type for single region in single ALLC

   .. py:method:: read_line(self, line)


   .. py:method:: summary(self)



.. py:function:: determine_datasets(regions, quantifiers, chrom_size_path, tmp_dir)


.. py:function:: count_single_region_set(allc_table, region_config, obs_dim, region_dim)

   Get cell-by-region-by-mc_types count matrix, save to zarr


.. py:function:: calculate_pv(data, reverse_value, obs_dim, var_dim, cutoff=0.9)


.. py:function:: count_single_zarr(allc_table, region_config, obs_dim, region_dim, output_path, obs_dim_dtype, count_dtype='uint32')

   process single region set and its quantifiers


.. py:function:: generate_dataset(allc_table, output_path, regions, quantifiers, chrom_size_path, obs_dim='cell', cpu=1, chunk_size=None)

   {generate_dataset_doc}

   :param allc_table: {allc_table_doc}
   :param output_path: Output path of the MCDS dataset
   :param regions: {regions_doc}
   :param quantifiers: {quantifiers_doc}
   :param chrom_size_path: {chrom_size_path_doc}
   :param obs_dim: {obs_dim_doc}
   :param cpu: {cpu_basic_doc}
   :param chunk_size: {chunk_size_doc}

   :returns:
   :rtype: output_path


