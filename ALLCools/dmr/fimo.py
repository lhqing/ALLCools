"""
Run FIMO to scan motif over FASTA sequences
"""
import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

from ..motif.utilities import split_meme_motif_file
from .get_fasta import get_fasta


def _fimo_runner(motif_path, fasta_path, output_path, path_to_fimo='',
                 raw_score_thresh=6, raw_p_value_thresh=1e-4, parse_genome_coords=False):
    """Run fimo for a single motif over single fasta file"""
    if path_to_fimo != '':
        path_to_fimo = path_to_fimo.rstrip('/') + '/'
    parse_genome_coords_opt = '--parse-genomic-coord' if parse_genome_coords else ''
    cmd = f'{path_to_fimo}fimo --text --skip-matched-sequence {parse_genome_coords_opt} ' \
          f'--thresh {raw_p_value_thresh} {motif_path} {fasta_path}'

    tmp_path = str(output_path) + "tmp"
    with open(tmp_path, 'w') as out:
        p = subprocess.Popen(shlex.split(cmd),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding='utf8')
        first = True
        while True:
            line = p.stdout.readline()
            if line == '':
                break
            if first:
                first = False
                out.write(line)
                continue
            score = line.split('\t')[-4]
            # when scan very small motifs, p-value is not very sig, but score can be different.
            if float(score) > raw_score_thresh:
                out.write(line)
        p.terminate()
        if p.returncode is not None and p.returncode != 0:
            raise OSError(f'{cmd} return code {p.returncode}')

    try:
        final_df = pd.read_csv(tmp_path, sep='\t')
        # same as sort -k1,1 -k2,2n, so in bedtools intersect, we can use --sorted
        if parse_genome_coords:
            final_df = final_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
        else:
            final_df = final_df.sort_values(by=['sequence_name', 'start']).reset_index(drop=True)
    except pd.errors.EmptyDataError:
        # if no result in tmp_path
        return ''

    # motif stats are not saved
    print(f'{motif_path}, N motif total={final_df.shape[0]}')
    final_df.to_hdf(output_path, key='data')
    subprocess.run(['rm', '-f', tmp_path], check=True)
    return output_path


def _scan_motif_over_fasta(meme_motif_file, fasta_path, output_dir, cpu=10, path_to_fimo='',
                           raw_score_thresh=6, raw_p_value_thresh=1e-4, parse_genome_coords=False):
    if isinstance(meme_motif_file, str):
        meme_motif_file = [meme_motif_file]

    output_dir = pathlib.Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    motif_file_records = split_meme_motif_file(meme_motif_file, output_dir)
    print(f'{len(motif_file_records)} motifs to count.')

    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for uid, name, motif_file_path in motif_file_records:
            future = executor.submit(_fimo_runner,
                                     path_to_fimo=path_to_fimo,
                                     motif_path=motif_file_path,
                                     fasta_path=fasta_path,
                                     output_path=output_dir / (pathlib.Path(motif_file_path).name[:-5] + '.bed.hdf'),
                                     raw_score_thresh=raw_score_thresh,
                                     raw_p_value_thresh=raw_p_value_thresh,
                                     parse_genome_coords=parse_genome_coords)
            futures[future] = (uid, name, motif_file_path)

        n = 0
        motif_records = []
        for future in as_completed(futures):
            output_path = future.result()
            uid, name, motif_file_path = futures[future]
            if output_path == '':
                # fimo return nothing, the parameter is too stringent for the motif.
                # by default parameters, this rarely happens.
                # Usually its the motif too small so p value is large.
                print(f'Motif {uid} {name} do not have any match under current settings.')
                continue
            motif_records.append((uid, name, str(motif_file_path), str(output_path)))
            n += 1
            if n % 100 == 0:
                print(f'Finish count {n} motifs.')

    lookup_table = pd.DataFrame(motif_records,
                                columns=['uid', 'name', 'motif_file_path', 'output_file_path']).set_index('uid')

    lookup_table.to_hdf(output_dir / f'LOOKUP_TABLE.hdf', key='data')
    return


def scan_motif_over_bed(bed_path, motif_file, genome_fasta,
                        output_dir, slop_b=None, chrom_size_path=None, use_region_name=False,
                        cpu=10, sort_mem_gbs=10, standard_length=500, sample_region=None,
                        raw_score_thresh=6, raw_p_value_thresh=1e-4):
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

    temp_fasta_path = str(output_dir / 'TEMP.fa')
    get_fasta([bed_path],
              genome_fasta,
              temp_fasta_path,
              slop_b=slop_b,
              chrom_size_path=chrom_size_path,
              use_region_name=use_region_name,
              cpu=cpu,
              sort_mem_gbs=sort_mem_gbs,
              standard_length=standard_length,
              merge=False,
              sample_region=sample_region,
              seed=1)

    _scan_motif_over_fasta(motif_file,
                           fasta_path=temp_fasta_path,
                           output_dir=output_dir, cpu=cpu, path_to_fimo='',
                           raw_score_thresh=raw_score_thresh,
                           raw_p_value_thresh=raw_p_value_thresh, parse_genome_coords=False)

    subprocess.run(['rm', '-f', temp_fasta_path])
    return
