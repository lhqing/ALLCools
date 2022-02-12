import logomaker
import numpy as np
import pandas as pd
import pathlib
import re
import shlex
import subprocess
from pybedtools import BedTool, cleanup

from .._open import open_gz
from ..utilities import parse_chrom_size


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


def meme_motif_file_to_dict(meme_motif_paths):
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
    return records


def single_meme_txt_to_pfm_df(text, bits_scale=True):
    sep = re.compile(r'[01].[0-9]+')
    enter = False
    pfm_rows = []
    for row in text.split('\n'):
        if enter:
            if row.startswith(' '):
                numbers = sep.findall(row)
                if len(numbers) == 4:
                    pfm_rows.append(list(map(float, numbers)))
        else:
            if row.startswith('letter-probability'):
                enter = True
            continue
    pfm = pd.DataFrame(pfm_rows, columns=['A', 'C', 'G', 'T'])

    if bits_scale:
        information_content = (pfm * np.log2(pfm / 0.25)).sum(axis=1)
        pfm = pfm.multiply(information_content, axis=0)
    return pfm


def meme_to_pfm_dict(meme_motif_paths, bits_scale=True):
    records = meme_motif_file_to_dict(meme_motif_paths)
    pfm_dict = {}
    for (uid, _), text in records.items():
        pfm_dict[uid] = single_meme_txt_to_pfm_df(text, bits_scale)
    return pfm_dict


def plot_pfm(pfm, ax=None, logo_kws=None):
    _logo_kws = dict(font_name='helvetica',
                     color_scheme=None,
                     vpad=.1,
                     ax=ax,
                     width=.8)
    if logo_kws is not None:
        _logo_kws.update(logo_kws)
    logo = logomaker.Logo(pfm, **_logo_kws)
    logo.ax.set_xlim([-1, len(pfm)])
    return logo


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
    records = meme_motif_file_to_dict(meme_motif_paths)

    motif_file_records = []
    for (uid, name), text in records.items():
        motif_file_path = output_dir / f'{uid}.meme'
        motif_file_records.append([uid, name, motif_file_path])
        with open(motif_file_path, 'w') as f:
            f.write(text)
    return motif_file_records
