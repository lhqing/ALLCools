from ALLCools.api import allc_to_region_count
from .utilities import *


def test_allc_to_region_count():
    for split_strand in [True, False]:
        for save_zero_cov in [True, False]:
            allc_to_region_count(allc_path=data_file_path('Cell01.allc.tsv.gz'),
                                 output_prefix=output_path(f'Cell01.allc.count.split{split_strand}.full{save_zero_cov}'),
                                 chrom_size_path=CHROM_SIZE_PATH,
                                 mc_contexts=['CGN', 'CHN'],
                                 split_strand=split_strand,
                                 region_bed_paths=[GENE_BED],
                                 region_bed_names=['gene'],
                                 bin_sizes=[1000, 100000],
                                 cov_cutoff=2,
                                 save_zero_cov=save_zero_cov,
                                 remove_tmp=False,
                                 cpu=1)
