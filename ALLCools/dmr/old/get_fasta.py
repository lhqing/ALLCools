import pathlib
import shlex
import subprocess

from pybedtools import BedTool, cleanup

from ALLCools._open import open_gz
from ALLCools.utilities import parse_chrom_size

import pandas as pd


def get_fasta(bed_file_paths, fasta_path, output_path, slop_b=None, chrom_size_path=None, use_region_name=False,
              cpu=1, sort_mem_gbs=1, standard_length=None, merge=False, sample_region=None, seed=1):
    """
    Extract genome sequence fasta using bed files

    Parameters
    ----------
    bed_file_paths
    fasta_path
    output_path
    slop_b
    chrom_size_path
    use_region_name
        If region names provided in the fourth column of bed file:
            if True: use region name as seq name
            else: use chr:start-end as seq name
    cpu
    sort_mem_gbs
    standard_length
    merge
    sample_region
    seed

    Returns
    -------

    """
    chrom_dict = None
    if chrom_size_path is not None:
        chrom_dict = parse_chrom_size(chrom_size_path)

    if isinstance(bed_file_paths, str):
        bed_file_paths = [bed_file_paths]
    output_path = str(pathlib.Path(output_path).resolve())

    temp_bed = output_path + '.tmp_input.bed'
    with open(temp_bed, 'w') as temp_f:
        for bed_file_path in bed_file_paths:
            if str(bed_file_path).endswith('gz'):
                opener = open_gz
            else:
                opener = open
            with opener(bed_file_path) as f:
                if standard_length is None:
                    temp_f.write(f.read())
                else:
                    if chrom_dict is None:
                        raise ValueError('chrom_size_path can not be None when standard_length is not None')
                    half = int(standard_length / 2)
                    for line in f:
                        ll = line.strip().split('\t')
                        center = (int(ll[1]) + int(ll[2])) / 2
                        ll[1] = str(int(max(center - half, 0)))
                        ll[2] = str(int(min(center + half, chrom_dict[ll[0]])))
                        temp_f.write('\t'.join(ll) + '\n')

    sorted_temp = output_path + '.tmp_sorted.bed'
    try:
        subprocess.run(
            shlex.split(
                f'sort -k1,1 -k2,2n --parallel={cpu} -S {sort_mem_gbs}G {temp_bed} -o {sorted_temp}'),
            stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)
    except subprocess.CalledProcessError:
        # old sort version don't have parallel option
        print('run sort without --parallel')
        try:
            subprocess.run(
                shlex.split(f'sort -k1,1 -k2,2n -S {sort_mem_gbs}G {temp_bed} -o {sorted_temp}'),
                stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            raise e

    sorted_bed = BedTool(sorted_temp)
    if slop_b:
        if chrom_size_path is None:
            raise ValueError('chrom_size_path can not be None when slop_b is not None')
        sorted_bed = sorted_bed.slop(b=slop_b, g=chrom_size_path)

    merged_temp = output_path + 'tmp_merge.bed'
    if merge:
        if use_region_name:
            print('can not use region name when merge is True')
            use_region_name = False
        sorted_bed.merge().moveto(merged_temp)
    else:
        sorted_bed.moveto(merged_temp)

    if sample_region is not None:
        bed_df = pd.read_csv(merged_temp, header=None, sep='\t')
        if sample_region <= bed_df.shape[0]:
            bed_df = bed_df.sample(sample_region, random_state=seed)
        bed_df.to_csv(merged_temp, sep='\t', index=None, header=None)

    name_option = '-name' if use_region_name else ''
    subprocess.run(
        shlex.split(f'bedtools getfasta -fi {fasta_path} -bed {merged_temp} -fo {output_path} {name_option}'),
        stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)

    subprocess.run(shlex.split(f'rm -f {temp_bed} {sorted_temp} {merged_temp}'))
    cleanup()
    return output_path
