import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from .utilities import get_fasta, split_meme_motif_file
from .._doc import *


def _single_ame_runner(fasta_path, motif_paths, ame_kws_str=''):
    if isinstance(motif_paths, list):
        motif_paths = ' '.join([str(i) for i in motif_paths])
    cmd = f'ame {ame_kws_str} {fasta_path} {motif_paths} '
    try:
        subprocess.run(shlex.split(cmd),
                       check=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       encoding='utf8')
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    return cmd


@doc_params(chrom_size_path=chrom_size_path_doc)
def ame(bed_file,
        motif_files,
        output_path,
        reference_fasta,
        background_files=None,
        cpu=1,
        standard_length=None,
        chrom_size_path=None,
        ame_kws=None):
    """\
    Motif enrichment analysis with AME from MEME Suite.

    Parameters
    ----------
    bed_file
        Single bed file input to calculate motif enrichment
    motif_files
        One or multiple motif files in MEME format
    output_path
        Output path of the AME result
    reference_fasta
        Reference fasta of the bed_file and background_files
    background_files
        One or multiple bed files contain region as the control set of motif enrichment
    cpu
        Number of CPUs to parallel AME
    standard_length
        If not None, will standard the input region and control region (if provided) in to same standard_length,
        each standardized region is centered to the original region center.
    chrom_size_path
        {chrom_size_path_doc}
        Must provide is standard_length is not None.
    ame_kws
        Additional options that will pass to AME, provide either a dict or string.
        If string provided, assume each option is separated by ';' and key-value is separated by '=',
        e.g. 'seed=1;method=fisher;scoring=avg'
    Returns
    -------

    """
    # test related software
    try:
        subprocess.run(['ame', '--version'], check=True)
    except subprocess.CalledProcessError as e:
        print('Make sure AME from MEME suite is correctly installed.')
        raise e
    try:
        subprocess.run(['bedtools', '--version'], check=True)  # used in get_fasta
    except subprocess.CalledProcessError as e:
        print('Make sure bedtools is correctly installed.')
        raise e

    output_path = str(output_path)

    output_dir = pathlib.Path(output_path + '_ame_temp')
    output_dir.mkdir(exist_ok=True)
    if isinstance(background_files, str):
        background_files = [background_files]

    # get bed FASTQ
    bed_file = pathlib.Path(bed_file)
    output_fasta_path = output_dir / (bed_file.name + '.fasta')
    get_fasta([bed_file],
              fasta_path=reference_fasta,
              output_path=output_fasta_path,
              standard_length=standard_length,
              slop_b=None,
              chrom_size_path=chrom_size_path,
              cpu=cpu, sort_mem_gbs=min(int(cpu * 3), 30))

    # get merged background FASTQ
    if background_files is None:
        control_fasta = '--shuffle--'
    else:
        control_fasta = output_dir / 'CONTROL_REGIONS.fasta'
        get_fasta(background_files,
                  fasta_path=reference_fasta,
                  output_path=control_fasta,
                  slop_b=None,
                  standard_length=standard_length,
                  chrom_size_path=chrom_size_path,
                  cpu=cpu, sort_mem_gbs=min(int(cpu * 3), 30))

    # split motif files based on cpu
    ame_motif_temp_dir = output_dir / 'MOTIF_TEMP'
    ame_motif_temp_dir.mkdir(exist_ok=True)
    motif_records = split_meme_motif_file(motif_files, ame_motif_temp_dir)
    motif_path_chunks = []
    step = len(motif_records) // cpu + 1
    for i in range(0, len(motif_records), step):
        # j is a list [motif_uid, motif_name, motif_path]
        # return by split_meme_motif_file
        motif_path_chunks.append([j[2] for j in motif_records[i:i + step]])

    # configure ame options
    ame_kws_str = ''
    _ame_kws = {'noseq': '',
                'verbose': 1,
                'control': control_fasta}
    if isinstance(ame_kws, dict):
        _ame_kws.update(ame_kws)
    elif isinstance(ame_kws, str):
        ame_kws_str += ame_kws
    else:
        pass
    for k, v in _ame_kws.items():
        ame_kws_str += f' --{k} {v}'

    # run single_bed_single_motif_runner
    with ProcessPoolExecutor(cpu) as pool:
        futures = []
        for chunk_id, chunk in enumerate(motif_path_chunks):
            _ame_kws['oc'] = output_dir / f'{chunk_id}_temp'
            future = pool.submit(_single_ame_runner,
                                 fasta_path=output_fasta_path,
                                 motif_paths=chunk,
                                 ame_kws_str=ame_kws_str)
            futures.append(future)

        for future in as_completed(futures):
            future.result()

    # merge results
    records = []
    for i in range(len(motif_path_chunks)):
        chunk_df = pd.read_csv(output_dir / f'{i}_temp/ame.tsv',
                               sep='\t', index_col=0, comment='#').reset_index(drop=True)
        records.append(chunk_df)
    total_records = pd.concat(records, sort=True).sort_values('E-value').reset_index(drop=True)
    total_records['motif_DB'] = total_records['motif_DB'].apply(lambda x: x.split('/')[-1])
    # add fold change
    total_records['TP/(FP+1)'] = total_records['TP'] / (total_records['FP'] + 1)
    total_records['log2(TP/(FP+1))'] = np.log2(total_records['TP/FP'])
    # save output
    total_records.to_csv(output_path, sep='\t')

    # clean temp dir
    subprocess.run(['rm', '-rf', output_dir], check=True)

    return total_records
