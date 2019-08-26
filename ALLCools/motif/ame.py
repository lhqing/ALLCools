import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

from .utilities import get_fasta, split_meme_motif_file


def _single_ame_runner(fasta_path, motif_paths, **ame_kws):
    options = ''
    for k, v in ame_kws.items():
        options += f' --{k} {v}'

    if isinstance(motif_paths, list):
        motif_paths = ' '.join([str(i) for i in motif_paths])
    cmd = f'ame {options} {fasta_path} {motif_paths} '
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


def parallel_ame(bed_file,
                 motif_files,
                 output_dir,
                 reference_fasta,
                 background_files=None,
                 cpu=1,
                 standard_length=None,
                 chrom_size_path=None,
                 ame_kws=None):
    """

    Run motif enrichment analysis on 1 bed_file with multiple motif_files,
    using ALL background_files (merged) as control region

    Parameters
    ----------
    bed_file
    motif_files
    output_dir
    reference_fasta
    background_files
    cpu
    standard_length
    chrom_size_path
    ame_kws

    Returns
    -------

    """
    output_dir = pathlib.Path(output_dir)
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

    # get merged background fastq
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
    _ame_kws = {'noseq': '',
                'verbose': 1,
                'control': control_fasta}
    if ame_kws is not None:
        _ame_kws.update(ame_kws)

    # run single_bed_single_motif_runner
    with ProcessPoolExecutor(cpu) as pool:
        futures = []
        for chunk_id, chunk in enumerate(motif_path_chunks):
            _ame_kws['oc'] = output_dir / f'{chunk_id}_temp'
            future = pool.submit(_single_ame_runner,
                                 fasta_path=output_fasta_path,
                                 motif_paths=chunk,
                                 **_ame_kws)
            futures.append(future)

        for future in as_completed(futures):
            future.result()

    # merge results
    records = []
    for i in range(len(motif_path_chunks)):
        chunk_df = pd.read_csv(output_dir / f'{i}_temp/ame.tsv',
                               sep='\t', index_col=0, comment='#').reset_index(drop=True)
        print(chunk_df.shape)
        records.append(chunk_df)
    total_records = pd.concat(records, sort=True)
    return total_records
