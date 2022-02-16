:py:mod:`ALLCools.mcds.mcds`
============================

.. py:module:: ALLCools.mcds.mcds


Module Contents
---------------

.. py:function:: make_obs_df_var_df(use_data, obs_dim, var_dim)


.. py:class:: MCDS(dataset, obs_dim=None, var_dim=None)

   Bases: :py:obj:`xarray.Dataset`

   MCDS Class

   .. py:attribute:: __slots__
      :annotation: = []

      

   .. py:method:: var_dim(self)
      :property:


   .. py:method:: obs_dim(self)
      :property:


   .. py:method:: obs_names(self)
      :property:


   .. py:method:: var_names(self)
      :property:


   .. py:method:: _verify_dim(self, dim, mode)


   .. py:method:: open(cls, mcds_paths, obs_dim='cell', use_obs=None, var_dim=None, chunks='auto', split_large_chunks=True, obj_to_str=True, engine=None)
      :classmethod:

      Take one or multiple MCDS file paths and create single MCDS concatenated on obs_dim

      :param mcds_paths: Single MCDS path or MCDS path pattern with wildcard or MCDS path list
      :param obs_dim: Dimension name of observations, default is 'cell'
      :param use_obs: Subset the MCDS by a list of observation IDs.
      :param var_dim: Which var_dim dataset to use, needed when MCDS has multiple var_dim stored in the same directory
      :param chunks: if not None, xarray will use chunks to load data as dask.array. The "auto" means xarray will
                     determine chunks automatically. For more options, read the `xarray.open_dataset` `chunks` parameter
                     documentation. If None, xarray will not use dask, which is not desired in most cases.
      :param split_large_chunks: Whether split large chunks in dask config array.slicing.split_large_chunks
      :param obj_to_str: Whether turn object coordinates into string data type
      :param engine: xarray engine used to store MCDS, if multiple MCDS provided, the engine need to be the same

      :returns:
      :rtype: MCDS


   .. py:method:: add_mc_frac(self, var_dim=None, da=None, normalize_per_cell=True, clip_norm_value=10, da_suffix='frac')

      Add posterior mC rate data array for certain feature type (var_dim).

      :param var_dim: Name of the feature type
      :param da: if None, will use f'{var_dim}_da'
      :param normalize_per_cell: if True, will normalize the mC rate data array per cell
      :param clip_norm_value: reset larger values in the normalized mC rate data array to this
      :param da_suffix: name suffix appended to the calculated mC rate data array


   .. py:method:: add_mc_rate(self, *args, **kwargs)


   .. py:method:: add_feature_cov_mean(self, obs_dim=None, var_dim=None, plot=True)

      Add feature cov mean across obs_dim.

      :param var_dim: Name of var dimension
      :param obs_dim: Name of obs dimension
      :param plot: If true, plot the distribution of feature cov mean

      :returns:
      :rtype: None


   .. py:method:: add_cell_metadata(self, metadata, obs_dim=None)


   .. py:method:: filter_feature_by_cov_mean(self, var_dim=None, min_cov=0, max_cov=999999)

      filter MCDS by feature cov mean. add_feature_cov_mean() must be called before this function.

      :param var_dim: Name of var dimension
      :param min_cov: Minimum cov cutoff
      :param max_cov: Maximum cov cutoff

      :returns:
      :rtype: MCDS


   .. py:method:: get_feature_bed(self, var_dim=None)

      Get a bed format data frame of the var_dim

      :param var_dim: Name of var_dim

      :returns:
      :rtype: pd.DataFrame


   .. py:method:: remove_black_list_region(self, black_list_path, var_dim=None, f=0.2)

      Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

      :param var_dim: Name of var_dim
      :param black_list_path: Path to the black list bed file
      :param f: Fraction of overlap when calling bedtools intersect

      :returns:
      :rtype: MCDS


   .. py:method:: remove_chromosome(self, exclude_chromosome, var_dim=None)

      Remove regions in specific chromosome

      :param var_dim: Name of var_dim
      :param exclude_chromosome: Chromosome to remove

      :returns:
      :rtype: MCDS (xr.Dataset)


   .. py:method:: calculate_hvf_svr(self, mc_type=None, var_dim=None, obs_dim=None, n_top_feature=5000, da_suffix='frac', plot=True)


   .. py:method:: calculate_hvf(self, mc_type=None, var_dim=None, obs_dim=None, min_disp=0.5, max_disp=None, min_mean=0, max_mean=5, n_top_feature=5000, bin_min_features=5, mean_binsize=0.05, cov_binsize=100, da_suffix='frac', plot=True)

      Calculate normalized dispersion to select highly variable features.

      :param mc_type: Type of mC to calculate
      :param var_dim: Name of variable
      :param obs_dim: Name of observation, default is cell
      :param min_disp: minimum dispersion for a feature to be considered
      :param max_disp: maximum dispersion for a feature to be considered
      :param min_mean: minimum mean for a feature to be considered
      :param max_mean: maximum mean for a feature to be considered
      :param n_top_feature: Top N feature to use as highly variable feature.
                            If set, all the cutoff will be ignored, HDF selected based on order of normalized dispersion.
      :param bin_min_features: Minimum number of features to be considered as a separate bin,
                               if bellow this number, the bin will be merged to its closest bin.
      :param mean_binsize: bin size to separate features across mean
      :param cov_binsize: bin size to separate features across coverage
      :param plot: If true, will plot mean, coverage and normalized dispersion scatter plots.

      :returns:
      :rtype: pd.DataFrame


   .. py:method:: get_score_adata(self, mc_type, quant_type, obs_dim=None, var_dim=None, sparse=True)


   .. py:method:: get_adata(self, mc_type=None, obs_dim=None, var_dim=None, da_suffix='frac', select_hvf=True, split_large_chunks=True)

      Get anndata from MCDS mC rate matrix
      :param mc_type: mC rate type
      :param var_dim: Name of variable
      :param da_suffix: Suffix of mC rate matrix
      :param obs_dim: Name of observation
      :param select_hvf: Select HVF or not, if True, will use mcds.coords['{var_dim}_{mc_type}_feature_select'] to select HVFs
      :param split_large_chunks: Whether split large chunks in dask config array.slicing.split_large_chunks

      :returns:
      :rtype: anndata.Anndata


   .. py:method:: merge_cluster(self, cluster_col, obs_dim=None, add_mc_frac=True, add_overall_mc=True, overall_mc_da='chrom100k_da')


   .. py:method:: to_region_ds(self, region_dim=None)


   .. py:method:: write_dataset(self, output_path, mode='w-', obs_dim=None, var_dims: Union[str, list] = None, use_obs=None, chunk_size=1000)

      Write MCDS into a on-disk zarr dataset. Data arrays for each var_dim will be saved in separate
      sub-directories of output_path.

      :param output_path: Path of the zarr dataset
      :param mode: 'w-' means write to output_path, fail if the path exists; 'w' means write to output_path,
                   overwrite if the var_dim sub-directory exists
      :param obs_dim: dimension name of observations
      :param var_dims: dimension name, or a list of dimension names of variables
      :param use_obs: Select AND order observations when write.
      :param chunk_size: The load and write chunks, set this as large as possible based on available memory.

      :returns:
      :rtype: output_path



