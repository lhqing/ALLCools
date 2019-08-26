import pathlib
import shlex
import subprocess

import numpy as np
from pybedtools import BedTool, cleanup
import pandas as pd
from .._open import open_gz
from ..utilities import parse_chrom_size


def get_fasta(bed_file_paths, fasta_path, output_path, slop_b=None, chrom_size_path=None,
              cpu=1, sort_mem_gbs=1, standard_length=None, merge=False, sample_region=None, seed=1):
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
            shlex.split(f'sort -k1,1 -k2,2n --parallel={cpu} -S {sort_mem_gbs}G {temp_bed} -o {sorted_temp}'),
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
        sorted_bed.merge().moveto(merged_temp)
    else:
        sorted_bed.moveto(merged_temp)

    if sample_region is not None:
        bed_df = pd.read_csv(merged_temp, header=None, sep='\t')
        if sample_region <= bed_df.shape[0]:
            bed_df = bed_df.sample(sample_region, random_state=seed)
        bed_df.to_csv(merged_temp, sep='\t', index=None, header=None)

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


def split_meme_motif_file(meme_motif_paths, output_dir):
    """
    Given multi motif meme format file, split into single motif meme format file

    Parameters
    ----------
    meme_motif_paths
    output_dir

    Returns
    -------

    """
    if isinstance(meme_motif_paths, str):
        meme_motif_paths = [meme_motif_paths]

    records = {}
    uid_set = set()
    for meme_motif_path in meme_motif_paths:
        with open(meme_motif_path) as f:
            first_line = f.readline()
            if not first_line.startswith('MEME version'):
                raise ValueError('Input file need to be MEME motif format.')

            f.seek(0)
            header = True
            header_text = ''
            motif_tmp_text = ''
            for line in f:
                if line.startswith('MOTIF'):
                    # save the previous one first
                    if motif_tmp_text != '':
                        records[(uid, name)] = header_text + motif_tmp_text
                    motif_tmp_text = line

                    try:
                        _, uid, name = line.strip('\n').split(' ')
                    except ValueError:
                        _, uid = line.strip('\n').split(' ')
                        name = ''
                    if uid in uid_set:
                        raise ValueError(f'Found duplicate motif uid {uid} in file {meme_motif_path}. '
                                         f'Motif uid should be unique across all meme files provided.')
                    else:
                        uid_set.add(uid)
                    header = False
                elif header:
                    header_text += line
                else:
                    motif_tmp_text += line
            if motif_tmp_text != '':
                records[(uid, name)] = header_text + motif_tmp_text

    motif_file_records = []
    for (uid, name), text in records.items():
        motif_file_path = output_dir / f'{uid}.meme'
        motif_file_records.append([uid, name, motif_file_path])
        with open(motif_file_path, 'w') as f:
            f.write(text)
    return motif_file_records
