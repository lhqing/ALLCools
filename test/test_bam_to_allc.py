import gzip

from ALLCools.api import bam_to_allc
from .utilities import *


def test_bam_to_allc_api():
    # default
    bam_to_allc(bam_path=data_file_path('Cell01.bam'),
                reference_fasta=REFERENCE_FASTA,
                output_path=output_path('Cell01.allc.tsv.gz'),
                cpu=1,
                num_upstr_bases=1,  # NOMe setting
                num_downstr_bases=2,
                min_mapq=10,
                min_base_quality=20,
                compress_level=5,
                save_count_df=True)

    with gzip.open(output_path('Cell01.allc.tsv.gz'), 'rt') as f:
        for line in f:
            ll = line.strip().split('\t')
            assert len(ll) == 7
            assert ll[2] in '+-'
            assert len(ll[3]) == 4
            assert ll[3][1] == 'C'
            assert int(ll[4]) <= int(ll[5])

    """
    # Parallel
    bam_to_allc(bam_path=data_file_path('Cell02.bam'),
                reference_fasta=REFERENCE_FASTA,
                output_path=output_path('Cell02.allc.tsv.gz'),
                cpu=5,
                num_upstr_bases=0,
                num_downstr_bases=2,
                min_mapq=10,
                min_base_quality=20,
                compress_level=5,
                save_count_df=False)

    with gzip.open(output_path('Cell02.allc.tsv.gz'), 'rt') as f:
        for line in f:
            ll = line.strip().split('\t')
            assert len(ll) == 7
            assert ll[2] in '+-'
            assert len(ll[3]) == 3
            assert ll[3][0] == 'C'
            assert int(ll[4]) <= int(ll[5])
    """
    return


def test_bam_to_allc_cli():
    cmd = f'allcools bam-to-allc ' \
          f'--bam_path {data_file_path("Cell03.bam")} ' \
          f'--reference_fasta {REFERENCE_FASTA} ' \
          f'--output_path {output_path("Cell03.allc.tsv.gz")} ' \
          f'--cpu 5 ' \
          f'--num_upstr_bases 0 ' \
          f'--num_downstr_bases 2 ' \
          f'--min_mapq 10 ' \
          f'--min_base_quality 20 ' \
          f'--compress_level 5 ' \
          f'--save_count_df'
    run_command(cmd)

    with gzip.open(output_path('Cell03.allc.tsv.gz'), 'rt') as f:
        for line in f:
            ll = line.strip().split('\t')
            assert len(ll) == 7
            assert ll[2] in '+-'
            assert len(ll[3]) == 3
            assert ll[3][0] == 'C'
            assert int(ll[4]) <= int(ll[5])

    # check if tabix works
    run_command(cmd=f"tabix {output_path('Cell03.allc.tsv.gz')} -l")
    return
