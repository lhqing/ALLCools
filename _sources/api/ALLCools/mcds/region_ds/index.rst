:py:mod:`ALLCools.mcds.region_ds`
=================================

.. py:module:: ALLCools.mcds.region_ds


Module Contents
---------------

.. py:function:: _bigwig_over_bed(bed: pandas.DataFrame, path, value_type='mean', dtype='float32')


.. py:function:: _region_bed_sorted(bed_path, g, bed_sorted)


.. py:function:: _bed_intersection(bed: pybedtools.BedTool, path, g, region_index, bed_sorted, fraction=0.2)


.. py:function:: _annotate_by_bigwigs_worker(dataset_path, region_dim, chrom_size_path, track_paths, output_path, dim, slop, value_type, dtype, **kwargs)


.. py:function:: _annotate_by_beds_worker(dataset_path, region_dim, chrom_size_path, slop, track_paths, dtype, dim, output_path, bed_sorted, fraction=0.2, **kwargs)


.. py:function:: _fisher_exact(row, alternative='two-sided')


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


   .. py:method:: from_bed(cls, bed, location, chrom_size_path, region_dim='region', sort_bed=True)
      :classmethod:

      Create empty RegionDS from a bed file.

      :param bed:
      :param location:
      :param region_dim:
      :param chrom_size_path:
      :param sort_bed:


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


   .. py:method:: annotate_by_beds(self, bed_table, dim, slop=100, chrom_size_path=None, chunk_size='auto', dtype='bool', bed_sorted=True, cpu=1, fraction=0.2, save=True)


   .. py:method:: get_feature(self, feature_name, dim=None, da_name=None)


   .. py:method:: scan_motifs(self, genome_fasta, cpu=1, standardize_length=500, motif_set_path=None, chrom_size_path=None, combine_cluster=True, fnr_fpr_fold=1000, chunk_size=None, dim='motif')


   .. py:method:: get_hypo_hyper_index(self, a, region_dim=None, region_state_da=None, sample_dim='sample', use_collapsed=True)


   .. py:method:: get_pairwise_differential_index(self, a, b, dmr_type='hypo', region_dim=None, region_state_da=None, sample_dim='sample', use_collapsed=True)


   .. py:method:: motif_enrichment(self, true_regions, background_regions, region_dim=None, motif_dim='motif-cluster', motif_da=None, alternative='two-sided')


   .. py:method:: sample_dmr_motif_enrichment(self, sample, region_dim=None, sample_dim='sample', motif_dim='motif-cluster', region_state_da=None, motif_da=None, alternative='two-sided', use_collapsed=True)


   .. py:method:: pairwise_dmr_motif_enrichment(self, a, b, dmr_type='hypo', region_dim=None, region_state_da=None, sample_dim='sample', motif_dim='motif-cluster', motif_da=None, alternative='two-sided')


   .. py:method:: object_coords_to_string(self, dtypes=None)


   .. py:method:: save(self, da_name=None, output_path=None, mode='w', change_region_dim=True)


   .. py:method:: get_coords(self, name)



