:py:mod:`ALLCools.mcds.utilities`
=================================

.. py:module:: ALLCools.mcds.utilities


Module Contents
---------------

.. py:data:: log
   

   

.. py:function:: calculate_posterior_mc_frac(mc_da, cov_da, var_dim=None, normalize_per_cell=True, clip_norm_value=10)


.. py:function:: calculate_posterior_mc_frac_lazy(mc_da, cov_da, var_dim, output_prefix, cell_chunk=20000, dask_cell_chunk=500, normalize_per_cell=True, clip_norm_value=10)

   Running calculate_posterior_mc_rate with dask array and directly save to disk.
   This is highly memory efficient. Use this for dataset larger then machine memory.

   :param mc_da:
   :param cov_da:
   :param var_dim:
   :param output_prefix:
   :param cell_chunk:
   :param dask_cell_chunk:
   :param normalize_per_cell:
   :param clip_norm_value:


.. py:function:: calculate_gch_rate(mcds, var_dim='chrom100k')


.. py:function:: get_mean_dispersion(x, obs_dim)


.. py:function:: highly_variable_methylation_feature(cell_by_feature_matrix, feature_mean_cov, obs_dim=None, var_dim=None, min_disp=0.5, max_disp=None, min_mean=0, max_mean=5, n_top_feature=None, bin_min_features=5, mean_binsize=0.05, cov_binsize=100)

   Adapted from Scanpy, the main difference is that,
   this function normalize dispersion based on both mean and cov bins.


.. py:function:: determine_engine(dataset_paths)


.. py:function:: obj_to_str(ds, coord_dtypes=None)


.. py:function:: write_ordered_chunks(chunks_to_write, final_path, append_dim, engine='zarr', coord_dtypes=None, dtype=None)


.. py:function:: convert_to_zarr(paths)

   Convert xarray.Dataset stored in other backends into zarr backend.


