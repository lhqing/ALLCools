import pathlib
import shlex
import subprocess

import numpy as np
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
    try:
        subprocess.run(
            shlex.split(f'sort -k1,1 -k2,2n --parallel={cpu} -S {sort_mem_gbs}G {temp_bed} -o {sorted_temp}'),
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
    sorted_bed.merge().moveto(merged_temp)

    subprocess.run(shlex.split(f'bedtools getfasta -fi {fasta_path} -bed {merged_temp} -fo {output_path}'),
                   stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)

    subprocess.run(shlex.split(f'rm -f {temp_bed} {sorted_temp} {merged_temp}'))
    cleanup()
    return output_path


def meme_to_homer(meme_path, homer_path, score_power=0.85):
    """
    Transfer MEME motif format into Homer motif format.
    Based on description here: http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html
    The score_power controls Log odds detection threshold by max_score ** score_power

    Parameters
    ----------
    meme_path
        Input path in meme format
    homer_path
        Output path in homer format
    score_power
        Log odds detection threshold = max_score ** score_power
    Returns
    -------

    """

    if score_power >= 1:
        raise ValueError('score_power must < 1')

    in_motif = False
    records = []
    with open(meme_path) as f:
        for line in f:
            if line.startswith('MOTIF'):
                if in_motif:
                    records.append(cur_char)
                    cur_char = ''
                else:
                    in_motif = True
                    cur_char = ''
            if not in_motif:
                continue
            else:
                cur_char += line

    with open(homer_path, 'w') as f:
        for record in records:
            lines = record.split('\n')
            head_line = lines[0]
            matrix_line = []
            enter_matrix = False
            for line in lines:
                if line.startswith('letter-probability matrix'):
                    enter_matrix = True
                    continue

                ll = line.strip().split('  ')
                ll = [i.strip() for i in ll]
                if enter_matrix and len(ll) == 4:
                    matrix_line.append(ll)
                else:
                    continue
            matrix = np.array(matrix_line).astype(float)
            best_score = np.log(matrix.max(axis=1) / 0.25).sum()
            uid = head_line.split(' ')[1]
            homer_head_line = f'>\t{uid}\t{best_score ** score_power}\n'
            matrix_line = '\n'.join(['\t'.join(line) for line in matrix_line]) + '\n'
            f.write(homer_head_line + matrix_line)
    return
