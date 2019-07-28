import pathlib
import shlex
import subprocess

from pybedtools import BedTool, cleanup

from .._open import open_gz


def get_fasta(bed_file_paths, fasta_path, output_path, slop_b=None, chrom_size_path=None,
              cpu=1, sort_mem_gbs=1):
    """
    Extract fasta using bed files
    The name of sequence in generated fasta file are: "chr:start-end"
    In fimo results, this name can be used to get original fasta position

    Parameters
    ----------
    bed_file_paths
    fasta_path
    output_path
    slop_b
    chrom_size_path
    cpu
    sort_mem_gbs

    Returns
    -------

    """
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
                temp_f.write(f.read())

    sorted_temp = output_path + '.tmp_sorted.bed'
    subprocess.run(shlex.split(f'sort -k1,1 -k2,2n --parallel={cpu} -S {sort_mem_gbs}G {temp_bed} -o {sorted_temp}'),
                   stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)

    sorted_bed = BedTool(sorted_temp)
    if slop_b:
        if chrom_size_path is None:
            raise ValueError('chrom_size_path can not be None when slop_b is not None')
        sorted_bed = sorted_bed.slop(b=slop_b, g=chrom_size_path)
    merged_temp = output_path + 'tmp_merge.bed'
    sorted_bed.merge().moveto(merged_temp)

    subprocess.run(shlex.split(f'bedtools getfasta -fi {fasta_path} -bed {merged_temp} -fo {output_path}'),
                   stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)

    subprocess.run(shlex.split(f'rm -f {temp_bed} {sorted_temp} {merged_temp}'))
    cleanup()
    return output_path
