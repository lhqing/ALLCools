from ALLCools.api import extract_allc
from .utilities import *


def test_extract_allc():
    for output_format in ['allc', 'bed5']:
        for strandness in ['both', 'split', 'merge']:
            for cpu in [1, 4]:
                extract_allc(allc_path=data_file_path('Cell01.allc.tsv.gz'),
                             output_prefix=output_path(f'Cell01.allc.extract.{output_format}_{strandness}_cpu{cpu}'),
                             mc_contexts=['CHN', 'CGN'],
                             chrom_size_path=CHROM_SIZE_PATH,
                             strandness=strandness,
                             output_format=output_format,
                             region=None,
                             cov_cutoff=9999,
                             tabix=True,
                             cpu=cpu)
    # TODO add correct answer and check the file carefully
