:py:mod:`ALLCools._allc_to_region_count`
========================================

.. py:module:: ALLCools._allc_to_region_count


Module Contents
---------------

.. py:function:: _bedtools_map(region_bed, site_bed, out_bed, chrom_size_path, save_zero_cov=True)

   Use bedtools map to map site_bed format into any bed file provided.


.. py:function:: _map_to_sparse_chrom_bin(site_bed, out_bed, chrom_size_path, bin_size=500)

   Calculate chromosome bins regional count, output is SPARSE,
   bin_id constructed from chrom_size_path and can be reproduce.


.. py:function:: allc_to_region_count(allc_path: str, output_prefix: str, chrom_size_path: str, mc_contexts: List[str], split_strand: bool = False, region_bed_paths: List[str] = None, region_bed_names: List[str] = None, bin_sizes: List[int] = None, cov_cutoff: int = 9999, save_zero_cov: bool = False, remove_tmp: bool = True, cpu: int = 1, binarize: bool = False)

   Calculate mC and cov at regional level. Region can be provided in 2 forms:
   1. BED file, provided by region_bed_paths, containing arbitrary regions and use bedtools map to calculate;
   2. Fix-size non-overlap genome bins, provided by bin_sizes.
   Form 2 is much faster to calculate than form 1.
   The output file is in 6-column bed-like format: chrom start end region_uid mc cov

   :param allc_path: {allc_path_doc}
   :param output_prefix: Path prefix of the output region count file.
   :param chrom_size_path: {chrom_size_path_doc}
   :param mc_contexts: {mc_contexts_doc}
   :param split_strand: {split_strand_doc}
   :param region_bed_paths: {region_bed_paths_doc}
   :param region_bed_names: {region_bed_names_doc}
   :param bin_sizes: {bin_sizes_doc}
   :param cov_cutoff: {cov_cutoff_doc}
   :param save_zero_cov: Whether to save the regions that have 0 cov, only apply to region count but not the chromosome count
   :param remove_tmp: Whether to remove the temporary BED file
   :param cpu: {cpu_basic_doc}
               This function parallel on region level at the extraction step
               and will generate a bunch of small files if cpu > 1.
               Do not use cpu > 1 for single cell region count.
               For single cell data, parallel on cell level is better.
   :param binarize: {binarize_doc}


.. py:function:: batch_allc_to_region_count(allc_series, output_dir, chrom_size_path, mc_contexts, split_strand, bin_sizes=None, region_bed_paths=None, region_bed_names=None, cov_cutoff=9999, cpu=5, binarize=False)


