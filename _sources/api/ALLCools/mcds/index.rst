:py:mod:`ALLCools.mcds`
=======================

.. py:module:: ALLCools.mcds

.. autoapi-nested-parse::

   Core Data Structure



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   correlation/index.rst
   mcds/index.rst
   region_ds/index.rst
   region_ds_utilities/index.rst
   utilities/index.rst


Package Contents
----------------

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



.. py:class:: RegionDS(dataset, region_dim=None, location=None, chrom_size_path=None)

   Bases: :py:obj:`xarray.Dataset`

   A multi-dimensional, in memory, array database.

   A dataset resembles an in-memory representation of a NetCDF file,
   and consists of variables, coordinates and attributes which
   together form a self describing dataset.

   Dataset implements the mapping interface with keys given by variable
   names and values given by DataArray objects for each variable name.

   One dimensional variables with name equal to their dimension are
   index coordinates used for label based indexing.

   To load data from a file or file-like object, use the `open_dataset`
   function.

   :param data_vars: A mapping from variable names to :py:class:`~xarray.DataArray`
                     objects, :py:class:`~xarray.Variable` objects or to tuples of
                     the form ``(dims, data[, attrs])`` which can be used as
                     arguments to create a new ``Variable``. Each dimension must
                     have the same length in all variables in which it appears.

                     The following notations are accepted:

                     - mapping {var name: DataArray}
                     - mapping {var name: Variable}
                     - mapping {var name: (dimension name, array-like)}
                     - mapping {var name: (tuple of dimension names, array-like)}
                     - mapping {dimension name: array-like}
                       (it will be automatically moved to coords, see below)

                     Each dimension must have the same length in all variables in
                     which it appears.
   :type data_vars: dict-like, optional
   :param coords: Another mapping in similar form as the `data_vars` argument,
                  except the each item is saved on the dataset as a "coordinate".
                  These variables have an associated meaning: they describe
                  constant/fixed/independent quantities, unlike the
                  varying/measured/dependent quantities that belong in
                  `variables`. Coordinates values may be given by 1-dimensional
                  arrays or scalars, in which case `dims` do not need to be
                  supplied: 1D arrays will be assumed to give index values along
                  the dimension with the same name.

                  The following notations are accepted:

                  - mapping {coord name: DataArray}
                  - mapping {coord name: Variable}
                  - mapping {coord name: (dimension name, array-like)}
                  - mapping {coord name: (tuple of dimension names, array-like)}
                  - mapping {dimension name: array-like}
                    (the dimension name is implicitly set to be the same as the
                    coord name)

                  The last notation implies that the coord name is the same as
                  the dimension name.
   :type coords: dict-like, optional
   :param attrs: Global attributes to save on this dataset.
   :type attrs: dict-like, optional

   .. rubric:: Examples

   Create data:

   >>> np.random.seed(0)
   >>> temperature = 15 + 8 * np.random.randn(2, 2, 3)
   >>> precipitation = 10 * np.random.rand(2, 2, 3)
   >>> lon = [[-99.83, -99.32], [-99.79, -99.23]]
   >>> lat = [[42.25, 42.21], [42.63, 42.59]]
   >>> time = pd.date_range("2014-09-06", periods=3)
   >>> reference_time = pd.Timestamp("2014-09-05")

   Initialize a dataset with multiple dimensions:

   >>> ds = xr.Dataset(
   ...     data_vars=dict(
   ...         temperature=(["x", "y", "time"], temperature),
   ...         precipitation=(["x", "y", "time"], precipitation),
   ...     ),
   ...     coords=dict(
   ...         lon=(["x", "y"], lon),
   ...         lat=(["x", "y"], lat),
   ...         time=time,
   ...         reference_time=reference_time,
   ...     ),
   ...     attrs=dict(description="Weather related data."),
   ... )
   >>> ds
   <xarray.Dataset>
   Dimensions:         (x: 2, y: 2, time: 3)
   Coordinates:
       lon             (x, y) float64 -99.83 -99.32 -99.79 -99.23
       lat             (x, y) float64 42.25 42.21 42.63 42.59
     * time            (time) datetime64[ns] 2014-09-06 2014-09-07 2014-09-08
       reference_time  datetime64[ns] 2014-09-05
   Dimensions without coordinates: x, y
   Data variables:
       temperature     (x, y, time) float64 29.11 18.2 22.83 ... 18.28 16.15 26.63
       precipitation   (x, y, time) float64 5.68 9.256 0.7104 ... 7.992 4.615 7.805
   .. attribute:: description

      Weather related data.

   Find out where the coldest temperature was and what values the
   other variables had:

   >>> ds.isel(ds.temperature.argmin(...))
   <xarray.Dataset>
   Dimensions:         ()
   Coordinates:
       lon             float64 -99.32
       lat             float64 42.21
       time            datetime64[ns] 2014-09-08
       reference_time  datetime64[ns] 2014-09-05
   Data variables:
       temperature     float64 7.182
       precipitation   float64 8.326
   .. attribute:: description

      Weather related data.

   .. py:attribute:: __slots__
      :annotation: = []

      

   .. py:method:: region_dim(self)
      :property:


   .. py:method:: chrom_size_path(self)
      :property:


   .. py:method:: location(self)
      :property:


   .. py:method:: from_bed(cls, bed, location, chrom_size_path, region_dim='region')
      :classmethod:

      Create empty RegionDS from a bed file.

      :param bed:
      :param location:
      :param region_dim:
      :param chrom_size_path:


   .. py:method:: open(cls, path, region_dim=None, use_regions=None, split_large_chunks=True, chrom_size_path=None, select_dir=None, engine='zarr')
      :classmethod:


   .. py:method:: _open_single_dataset(cls, path, region_dim, split_large_chunks=True, chrom_size_path=None, location=None, engine=None)
      :classmethod:

      Take one or multiple RegionDS file paths and create single RegionDS concatenated on region_dim

      :param path: Single RegionDS path or RegionDS path pattern with wildcard or RegionDS path list
      :param region_dim: Dimension name of regions
      :param split_large_chunks: Split large dask array chunks if true

      :returns:
      :rtype: RegionDS


   .. py:method:: iter_index(self, chunk_size=100000, dim=None)


   .. py:method:: iter_array(self, chunk_size=100000, dim=None, da=None, load=False)


   .. py:method:: get_fasta(self, genome_fasta, output_path, slop=None, chrom_size_path=None, standardize_length=None)


   .. py:method:: get_bed(self, with_id=True, bedtools=False, slop=None, chrom_size_path=None, standardize_length=None)


   .. py:method:: _chunk_annotation_executor(self, annotation_function, cpu, save=True, **kwargs)


   .. py:method:: annotate_by_bigwigs(self, bigwig_table, dim, slop=100, chrom_size_path=None, value_type='mean', chunk_size='auto', dtype='float32', cpu=1, save=True)


   .. py:method:: annotate_by_bed(self, bed_table, dim, slop=100, chrom_size_path=None, chunk_size='auto', dtype='bool', bed_sorted=True, cpu=1, save=True)


   .. py:method:: get_feature(self, feature_name, dim=None, da_name=None)


   .. py:method:: scan_motifs(self, genome_fasta, cpu=1, standardize_length=500, motif_set_path=None, chrom_size_path=None, combine_cluster=True, fnr_fpr_fold=1000, chunk_size=None, dim='motif')


   .. py:method:: get_hypo_hyper_index(self, a, region_dim=None, region_state_da=None, sample_dim='sample', use_collapsed=True)


   .. py:method:: get_pairwise_differential_index(self, a, b, dmr_type='hypo', region_dim=None, region_state_da=None, sample_dim='sample', use_collapsed=True)


   .. py:method:: motif_enrichment(self, true_regions, background_regions, region_dim=None, motif_dim='motif-cluster', motif_da=None, alternative='two-sided')


   .. py:method:: sample_dmr_motif_enrichment(self, sample, region_dim=None, sample_dim='sample', motif_dim='motif-cluster', region_state_da=None, motif_da=None, alternative='two-sided', use_collapsed=True)


   .. py:method:: pairwise_dmr_motif_enrichment(self, a, b, dmr_type='hypo', region_dim=None, region_state_da=None, sample_dim='sample', motif_dim='motif-cluster', motif_da=None, alternative='two-sided')


   .. py:method:: object_coords_to_string(self, dtypes=None)


   .. py:method:: save(self, output_path=None, mode='w', change_region_dim=True)


   .. py:method:: get_coords(self, name)



