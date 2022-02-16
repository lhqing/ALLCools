:py:mod:`ALLCools.sandbox.motif.motif_scan`
===========================================

.. py:module:: ALLCools.sandbox.motif.motif_scan


Module Contents
---------------

.. py:function:: _count_data_to_xarray(total_data: dict, allc_series)

   total_data: key is ALLC path, value is series of motif cov/mc count sum from all regions of that ALLC


.. py:function:: _count_single_cmotif_bin_on_multiple_allc(cmotif_dict_path, allc_paths, region, count_binary, context_to_pattern)


.. py:function:: allc_motif_scan(allc_table, output_path, mc_contexts, c_motif_dir, count_binary=True, cpu=1)

   Scan a list of ALLC files using a C-Motif database.
   C-Motif Database, can be generated via 'allcools generate-cmotif-database'
   Save the integrated multi-dimensional array into netCDF4 format using xarray.

   :param allc_table: {allc_table_doc}
   :param mc_contexts:
   :param c_motif_dir: A directory contains list of .msg files, each file records a dict,
                       position is key, value is tuple of motif direction and id
   :param output_path:
   :param count_binary: Only use this for single cell allc, instead of sum mC or cov directly,
                        will transfer mC and cov into [0, 1] when there is not ambiguity.
   :param cpu: {cpu_basic_doc}


