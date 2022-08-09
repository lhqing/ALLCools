:py:mod:`ALLCools.motif.snakemake`
==================================

.. py:module:: ALLCools.motif.snakemake


Module Contents
---------------

.. py:function:: prepare_motif_scan_snakemake(output_dir, fasta_path, region_dim, motif_dim='motif', motif_set_path=None, chunk_size=100000, combine_cluster=True, fnr_fpr_fold=1000, cpu=10)

   Prepare snakemake rules for motif scan.


.. py:function:: check_snakemake_success(output_dir)

   Check if snakemake jobs are finished.


.. py:function:: save_motif_chunks(motif_chunk_dir, region_dim, output_path, is_motif_cluster)

   Save motif chunks to xarray/zarr.


