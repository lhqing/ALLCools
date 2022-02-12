import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from .utilities import split_meme_motif_file
from ..dmr.get_fasta import get_fasta
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


@doc_params(chrom_size_path_doc=chrom_size_path_doc)
def ame(bed_file,
        motif_files,
        output_path,
        reference_fasta,
        background_files=None,
        cpu=1,
        standard_length=None,
        chrom_size_path=None,
        ame_kws=None,
        sample=None,
        seed=1):
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
        If string provided, will directly insert into the final AME command,
        see AME documentation about AME options.
        e.g. '--seed 1 --method fisher --scoring avg'
    sample
        Sample input regions to this number to reduce run time or test
    seed
        Seed for random sample input regions, only apply when sample is not None
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
    df = pd.read_csv(bed_file, sep='\t', header=None)
    print(df.shape[0], 'total regions')
    if sample is not None and sample <= df.shape[0]:
        df = df.sample(sample, random_state=seed)
        df.to_csv(output_dir / 'INPUT_BED_SAMPLED.bed',
                  sep='\t', index=None, header=None)
        bed_file = output_dir / 'INPUT_BED_SAMPLED.bed'
    n_region = df.shape[0]
    print(n_region, 'input regions')

    bed_file = pathlib.Path(bed_file)
    output_fasta_path = output_dir / (bed_file.name + '.fasta')
    get_fasta([bed_file],
              fasta_path=reference_fasta,
              output_path=output_fasta_path,
              standard_length=standard_length,
              slop_b=None,
              chrom_size_path=chrom_size_path,
              cpu=cpu,
              merge=False,
              sort_mem_gbs=min(int(cpu * 3), 30))

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
                  cpu=cpu,
                  merge=False,
                  sample_region=min(100000, n_region),
                  sort_mem_gbs=min(int(cpu * 3), 30))

    # split motif files based on cpu
    ame_motif_temp_dir = output_dir / 'MOTIF_TEMP'
    ame_motif_temp_dir.mkdir(exist_ok=True)
    motif_records = split_meme_motif_file(motif_files, ame_motif_temp_dir)
    motif_path_chunks = []
    step = min(len(motif_records) // cpu + 1, 10)
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
            _ame_kws_str = ame_kws_str
            chunk_output_dir = output_dir / f'{chunk_id}_temp'
            _ame_kws_str += f' --oc {chunk_output_dir}'
            future = pool.submit(_single_ame_runner,
                                 fasta_path=output_fasta_path,
                                 motif_paths=chunk,
                                 ame_kws_str=_ame_kws_str)
            futures.append(future)

        for future in as_completed(futures):
            future.result()

    # merge results
    records = []
    for i in range(len(motif_path_chunks)):
        ame_out_path = output_dir / f'{i}_temp/ame.tsv'
        # some rare cases AME output file has no header,I have no idea why
        # seems to be an AME bug.
        # So here I read the first line to check if the file has header
        with open(ame_out_path) as f:
            has_header = f.readline().startswith('rank')

        try:
            if has_header:
                chunk_df = pd.read_csv(ame_out_path,
                                       sep='\t', index_col=0, comment='#').reset_index(drop=True)
            else:
                chunk_df = pd.read_csv(ame_out_path, header=None,
                                       sep='\t', index_col=0, comment='#').reset_index(drop=True)
                chunk_df.columns = ['motif_DB', 'motif_ID', 'motif_alt_ID', 'consensus', 'p-value',
                                    'adj_p-value', 'E-value', 'tests', 'FASTA_max', 'pos', 'neg', 'PWM_min',
                                    'TP', '%TP', 'FP', '%FP']
                chunk_df.index.name = 'rank'

            print(i, chunk_df.shape)
        except pd.errors.EmptyDataError:
            continue
        records.append(chunk_df)
    total_records = pd.concat(records, sort=True)
    total_records = total_records[total_records['motif_DB'].notna()]
    total_records = total_records[total_records['%TP'].notna()]

    total_records['E-value'] = total_records['E-value'].fillna(100).astype(float)
    total_records = total_records.sort_values('E-value').reset_index(drop=True)
    total_records['motif_DB'] = total_records['motif_DB'].apply(
        lambda x: x.split('/')[-1] if isinstance(x, str) else x)

    # add fold change
    total_records['TP/FP'] = (total_records['%TP'] + 0.0001) / (total_records['%FP'] + 0.0001)
    total_records['log2(TP/FP)'] = np.log2(total_records['TP/FP'])
    # save output
    total_records.to_csv(output_path, sep='\t')

    # clean temp dir
    subprocess.run(['rm', '-rf', output_dir], check=True)

    return total_records
