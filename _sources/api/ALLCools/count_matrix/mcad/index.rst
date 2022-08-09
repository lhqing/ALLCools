:py:mod:`ALLCools.count_matrix.mcad`
====================================

.. py:module:: ALLCools.count_matrix.mcad


Module Contents
---------------

.. py:function:: _read_region_bed(bed_path)


.. py:function:: bin_sf(cov, mc, p)


.. py:function:: _count_single_allc(allc_path, bed_path, mc_pattern, output_dir, cutoff=0.9, reverse_value=False)


.. py:function:: generate_mcad(allc_table, bed_path, output_prefix, mc_context, cpu=1, cleanup=True, cutoff=0.9, reverse_value=False)

   Generate MCAD from ALLC files.

   :param allc_table: {allc_table_doc}
   :param bed_path: {bed_path_doc}
   :param cpu: {cpu_doc}
   :param output_prefix: Output prefix of the MCAD, a suffix ".mcad" will be added.
   :param mc_context: {mc_context_doc}
   :param cleanup: Whether remove temp files or not
   :param cutoff: Values smaller than cutoff will be stored as 0, which reduces the file size
   :param reverse_value: If true, use cdf instead of sf to make hyper-methylation events having higher values


