:py:mod:`ALLCools.mcds.mcds`
============================

.. py:module:: ALLCools.mcds.mcds


Module Contents
---------------

.. py:class:: MCDS(dataset)

   Bases: :py:obj:`xarray.Dataset`

   MCDS Class

   .. py:attribute:: __slots__
      :annotation: = []

      

   .. py:method:: open(cls, mcds_paths, obs_dim='cell', use_obs=None, split_large_chunks=True)
      :classmethod:

      Take one or multiple MCDS file paths and create single MCDS concatenated on obs_dim

      :param mcds_paths: Single MCDS path or MCDS path pattern with wildcard or MCDS path list
      :param obs_dim: Dimension name of observations, default is 'cell'
      :param use_obs: Subset the MCDS by a list of observation IDs.
      :param split_large_chunks: Whether split large chunks in dask config array.slicing.split_large_chunks

      :returns:
      :rtype: MCDS


   .. py:method:: add_mc_frac(self, var_dim, da=None, normalize_per_cell=True, clip_norm_value=10, da_suffix='frac')

      Add posterior mC rate data array for certain feature type (var_dim).

      :param var_dim: Name of the feature type
      :param da: if None, will use f'{var_dim}_da'
      :param normalize_per_cell: if True, will normalize the mC rate data array per cell
      :param clip_norm_value: reset larger values in the normalized mC rate data array to this
      :param da_suffix: name suffix appended to the calculated mC rate data array


   .. py:method:: add_mc_rate(self, *args, **kwargs)


   .. py:method:: add_feature_cov_mean(self, var_dim, obs_dim='cell', plot=True)

      Add feature cov mean across obs_dim.

      :param var_dim: Name of var dimension
      :param obs_dim: Name of obs dimension
      :param plot: If true, plot the distribution of feature cov mean

      :returns:
      :rtype: None


   .. py:method:: add_cell_metadata(self, metadata, obs_name='cell')


   .. py:method:: filter_feature_by_cov_mean(self, var_dim, min_cov=0, max_cov=999999)

      filter MCDS by feature cov mean. add_feature_cov_mean() must be called before this function.

      :param var_dim: Name of var dimension
      :param min_cov: Minimum cov cutoff
      :param max_cov: Maximum cov cutoff

      :returns:
      :rtype: MCDS


   .. py:method:: get_feature_bed(self, var_dim)

      Get a bed format data frame of the var_dim

      :param var_dim: Name of var_dim

      :returns:
      :rtype: pd.DataFrame


   .. py:method:: remove_black_list_region(self, var_dim, black_list_path, f=0.2)

      Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path

      :param var_dim: Name of var_dim
      :param black_list_path: Path to the black list bed file
      :param f: Fraction of overlap when calling bedtools intersect

      :returns:
      :rtype: MCDS


   .. py:method:: remove_chromosome(self, var_dim, exclude_chromosome)

      Remove regions in specific chromosome

      :param var_dim: Name of var_dim
      :param exclude_chromosome: Chromosome to remove

      :returns:
      :rtype: MCDS (xr.Dataset)


   .. py:method:: calculate_hvf_svr(self, var_dim, mc_type, obs_dim='cell', n_top_feature=5000, da_suffix='frac', plot=True)


   .. py:method:: calculate_hvf(self, mc_type, var_dim, obs_dim='cell', min_disp=0.5, max_disp=None, min_mean=0, max_mean=5, n_top_feature=5000, bin_min_features=5, mean_binsize=0.05, cov_binsize=100, plot=True)

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


   .. py:method:: get_adata(self, mc_type, var_dim, da_suffix='frac', obs_dim='cell', select_hvf=True, split_large_chunks=True)

      Get anndata from MCDS mC rate matrix
      :param mc_type: mC rate type
      :param var_dim: Name of variable
      :param da_suffix: Suffix of mC rate matrix
      :param obs_dim: Name of observation
      :param select_hvf: Select HVF or not, if True, will use mcds.coords['{var_dim}_{mc_type}_feature_select'] to select HVFs
      :param split_large_chunks: Whether split large chunks in dask config array.slicing.split_large_chunks

      :returns:
      :rtype: anndata.Anndata


   .. py:method:: merge_cluster(self, cluster_col, obs_dim='cell', add_mc_frac=True, add_overall_mc=True, overall_mc_da='chrom100k_da')


   .. py:method:: to_region_ds(self, region_dim)



