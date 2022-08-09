:py:mod:`ALLCools.mcds.mcds`
============================

.. py:module:: ALLCools.mcds.mcds


Module Contents
---------------

.. py:function:: _make_obs_df_var_df(use_data, obs_dim, var_dim)


.. py:class:: MCDS(dataset, coords=None, attrs=None, obs_dim=None, var_dim=None)

   Bases: :py:obj:`xarray.Dataset`

   The MCDS Class.

   .. py:attribute:: __slots__
      :annotation: = []

      

   .. py:method:: var_dim()
      :property:

      Name of the feature dimension.


   .. py:method:: obs_dim()
      :property:

      Name of the observation dimension.


   .. py:method:: obs_names()
      :property:

      Get obs_names.


   .. py:method:: var_names()
      :property:

      Get var_names.


   .. py:method:: _verify_dim(dim, mode)


   .. py:method:: get_var_dims(mcds_paths)
      :classmethod:

      Get var_dim from MCDS files.


   .. py:method:: open(mcds_paths, obs_dim='cell', use_obs=None, var_dim=None, chunks='auto', split_large_chunks=False, obj_to_str=True, engine=None, coords='minimal', compat='override', **kwargs)
      :classmethod:

      Take one or multiple MCDS file paths and create single MCDS concatenated on obs_dim.

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
      :param coords: the coords parameter of :py:func:`xarray.open_mfdataset` function,
                     default is "minimal", means only coordinates in which the dimension already appears are included.
      :param compat: the compat parameter of :py:func:`xarray.open_mfdataset` function,
                     default is "override", means skip comparing variables with the same name and pick variable from first MCDS.
      :param kwargs: Additional arguments passed on to :py:func:`xarray.open_dataset`
                     or :py:func:`xarray.open_mfdataset` function.

      :rtype: MCDS


   .. py:method:: add_mc_frac(var_dim=None, da=None, normalize_per_cell=True, clip_norm_value=10, da_suffix='frac')

      Add posterior mC rate data array for certain feature type (var_dim).

      :param var_dim: Name of the feature type
      :param da: if None, will use f'{var_dim}_da'
      :param normalize_per_cell: if True, will normalize the mC rate data array per cell
      :param clip_norm_value: reset larger values in the normalized mC rate data array to this
      :param da_suffix: name suffix appended to the calculated mC rate data array


   .. py:method:: _calculate_frac(var_dim, da, normalize_per_cell, clip_norm_value)

      Calculate mC frac data array for certain feature type (var_dim).


   .. py:method:: add_m_value(var_dim=None, da=None, alpha=0.01, normalize_per_cell=True, da_suffix='mvalue')

      Add m value data array for certain feature type (var_dim).

      M-Value is a transformation of the posterior mC fraction data array to a log ratio scale.
      M = np.log2((frac + alpha) / (1 - frac + alpha)).

      :param var_dim: Name of the feature type
      :param da: DataArray name. if None, will use f'{var_dim}_da'
      :param alpha: alpha value for the transformation regularization
      :param normalize_per_cell: if True, will normalize the mC rate data array per cell
      :param da_suffix: name suffix appended to the calculated mC rate data array


   .. py:method:: add_mc_rate(*args, **kwargs)

      Add mC fraction data array (Deprecated).


   .. py:method:: add_feature_cov_mean(obs_dim=None, var_dim=None, plot=True, da_name=None)

      Add feature cov mean across obs_dim.

      :param var_dim: Name of var dimension
      :param obs_dim: Name of obs dimension
      :param plot: If true, plot the distribution of feature cov mean
      :param da_name: Name of the calculated data array, if None, will use f'{var_dim}_da'


   .. py:method:: add_cell_metadata(metadata, obs_dim=None)

      Add cell metadata table to the MCDS.


   .. py:method:: filter_feature_by_cov_mean(var_dim=None, min_cov=0, max_cov=999999)

      Filter MCDS by feature cov mean. add_feature_cov_mean() must be called before this function.

      :param var_dim: Name of var dimension
      :param min_cov: Minimum cov cutoff
      :param max_cov: Maximum cov cutoff

      :rtype: MCDS


   .. py:method:: get_feature_bed(var_dim=None)

      Get a bed format data frame of the var_dim.

      :param var_dim: Name of var_dim

      :rtype: pd.DataFrame


   .. py:method:: remove_black_list_region(black_list_path, var_dim=None, f=0.2)

      Remove regions overlap (bedtools intersect -f {f}) with regions in the black_list_path.

      :param var_dim: Name of var_dim
      :param black_list_path: Path to the black list bed file
      :param f: Fraction of overlap when calling bedtools intersect

      :rtype: MCDS


   .. py:method:: remove_chromosome(exclude_chromosome=None, include_chromosome=None, var_dim=None)

      Remove regions in specific chromosome.

      :param var_dim: Name of var_dim
      :param exclude_chromosome: if provided, only these chromosomes will be removed
      :param include_chromosome: if provided, only these chromosomes will be kept

      :rtype: MCDS (xr.Dataset)


   .. py:method:: _get_da_name(var_dim, da_suffix)


   .. py:method:: calculate_hvf_svr(mc_type=None, var_dim=None, obs_dim=None, n_top_feature=5000, da_name=None, da_suffix='frac', plot=True)

      Calculate the highly variable features (hvf) with the Support Vector Regression model.


   .. py:method:: calculate_hvf(mc_type=None, var_dim=None, obs_dim=None, min_disp=0.5, max_disp=None, min_mean=0, max_mean=5, n_top_feature=5000, bin_min_features=5, mean_binsize=0.05, cov_binsize=100, da_name=None, da_suffix='frac', plot=True)

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
      :param da_name: Name of da to use, default is None, infer from var_dim and da_suffix
      :param da_suffix: Suffix to add to the name of the dataarray
      :param plot: If true, will plot mean, coverage and normalized dispersion scatter plots.

      :rtype: pd.DataFrame


   .. py:method:: get_count_adata(da_name, obs_dim=None, var_dim=None, sparse=True, loading_chunk=10000, binarize_cutoff=None, dtype='float32', use_vars=None, split_large_chunks=False)

      Convert a cell-by-feature count dataarray to an adata object.

      :param da_name: Name of the dataarray
      :param obs_dim: Name of observation
      :param var_dim: Name of variable
      :param sparse: If true, will adata.X convert to sparse matrix
      :param loading_chunk: Chunk size to load data in memory
      :param binarize_cutoff: If not None, will binarize the dataarray with values > binarize_cutoff as 1,
                              values <= binarize_cutoff as 0
      :param dtype: The final dtype of the adata.X, if None, will use the dtype of the dataarray
      :param use_vars: If not None, will use the specified variables as vars
      :param split_large_chunks: dask array.slicing.split_large_chunks parameter

      :rtype: anndata.AnnData


   .. py:method:: get_score_adata(mc_type, quant_type, obs_dim=None, var_dim=None, sparse=True, dtype='float32', loading_chunk=50000, binarize_cutoff=None)

      Convert a cell-by-feature methylation score dataarray to an adata object.

      :param mc_type: Name of the methylation type
      :param quant_type: Name of the quantification type, can be "hypo-score" or "hyper-score"
      :param obs_dim: Name of observation
      :param var_dim: Name of variable
      :param sparse: If true, will convert adata.X to sparse matrix
      :param dtype: If not None, will use the dtype of the adata.X
      :param loading_chunk: Chunk size to load data in memory
      :param binarize_cutoff: If not None, will binarize the dataarray with values > binarize_cutoff as 1,
                              values <= binarize_cutoff as 0

      :rtype: anndata.AnnData


   .. py:method:: add_feature_selection_column(feature_select, col_name='VAR_DIM_feature_select', var_dim=None)

      Manually add a feature selection column to the MCDS.


   .. py:method:: get_adata(mc_type=None, obs_dim=None, var_dim=None, da_name=None, da_suffix='frac', select_hvf=True, dtype='float32', split_large_chunks=False)

      Get anndata from MCDS mC rate matrix.

      :param mc_type: mC rate type
      :param var_dim: Name of variable
      :param da_suffix: Suffix of mC rate matrix
      :param obs_dim: Name of observation
      :param select_hvf: Select HVF or not, if True, will use mcds.coords['{var_dim}_{mc_type}_feature_select'] to select HVFs
      :param dtype: data type of adata.X
      :param split_large_chunks: Whether split large chunks in dask config array.slicing.split_large_chunks

      :rtype: anndata.Anndata


   .. py:method:: merge_cluster(cluster_col, obs_dim=None, add_mc_frac=True, add_overall_mc=True, overall_mc_da='chrom100k_da')

      Merge cell MCDS into cluster MCDS by sum on the obs_dim.


   .. py:method:: to_region_ds(region_dim=None)

      Turn the MCDS into a RegionDS.


   .. py:method:: write_dataset(output_path, mode='w-', obs_dim=None, var_dims: Union[str, list] = None, use_obs=None, chunks='auto')

      Write MCDS into an on-disk zarr dataset.

      Data arrays for each var_dim will be saved in separate sub-directories of output_path.
      The use_obs can be used to select and order observation accordingly.

      :param output_path: Path of the zarr dataset
      :param mode: 'w-' means write to output_path, fail if the path exists; 'w' means write to output_path,
                   overwrite if the var_dim sub-directory exists
      :param obs_dim: dimension name of observations
      :param var_dims: dimension name, or a list of dimension names of variables
      :param use_obs: Select and order observations according to this parameter when wrote to output_path.
      :param chunks: Zarr chunks on disk and in memory.

      :rtype: output_path


   .. py:method:: _write_dataset(output_path, mode='w-', obs_dim=None, var_dims: Union[str, list] = None, use_obs=None, chunks='auto')


   .. py:method:: save_feature_chunk_data(da_name, output_zarr_path, da_suffix='_fc', var_dim=None, loading_chunk=1000, var_chunk_size=1, obs_chunk_size=500000, compress_level=1, dtype=None)

      Save a data array to zarr dataset, which is chunked along the var_dim.

      This zarr dataset is useful when loading data from one or several specific
      features, such as making a gene plot.

      :param da_name: Name of data array to save.
      :param output_zarr_path: Path to output zarr dataset.
      :param da_suffix: Suffix to add to the name of the data array.
      :param var_dim: Name of var_dim. If None, use self.var_dim.
      :param loading_chunk: Number of var to load at a time.
      :param var_chunk_size: the var_dim chunk size of the output zarr dataset.
      :param obs_chunk_size: the obs_dim chunk size of the output zarr dataset.
      :param compress_level: the compress level of the output zarr dataset.
      :param dtype: the dtype of the output zarr dataset.



