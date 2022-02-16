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

   Generate multiple methylation datasets with a set of allc_table,
   a list of region sets and quantifiers for each region set.

   :param allc_table:
   :param output_path:
   :param regions:
   :param quantifiers:
   :param chrom_size_path:
   :param obs_dim:
   :param cpu:
   :param chunk_size:

   :returns:
   :rtype: output_path


