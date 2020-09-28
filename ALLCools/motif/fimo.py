"""
Run FIMO to scan motif over FASTA sequences
"""
import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

from .utilities import split_meme_motif_file

"""
Run FIMO to scan motif over FASTA sequences
"""


def _fimo_runner(motif_path, fasta_path, output_path, path_to_fimo='',
                 raw_score_thresh=6, raw_p_value_thresh=5e-4, is_genome_fasta=False,
                 top_n=300000):
    """Run fimo for a single motif over single fasta file"""
    if path_to_fimo != '':
        path_to_fimo = path_to_fimo.rstrip('/') + '/'
    cmd = f'{path_to_fimo}fimo --text --skip-matched-sequence --thresh {raw_p_value_thresh} {motif_path} {fasta_path}'

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
            score = line.split('\t')[6]
            # when scan very small motifs, p-value is not very sig, but score can be different.
            if float(score) > raw_score_thresh:
                out.write(line)
        p.terminate()
        if p.returncode is not None and p.returncode != 0:
            raise OSError(f'{cmd} returncode {p.returncode}')

    try:
        fimo_out = pd.read_csv(tmp_path, sep='\t')
        # same as sort -k1,1 -k2,2n, so in bedtools intersect, we can use --sorted
        fimo_out = fimo_out.sort_values(by=['sequence_name', 'start']).reset_index(drop=True)
        """
        # p-value correction and filtering
        # this is not very useful when scan large scale regions.
        
        reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(
            fimo_out['p-value'],
            alpha=alpha,
            method=correct_method,
            is_sorted=False,
            returnsorted=False)
        print(f'{motif_path}, alpha={alpha}, N motif total={reject.size}, '
              f'N motif remain={reject.sum()}, alphacBonf={alphacBonf}')
        fimo_out = fimo_out[reject].copy()
        """

    except pd.errors.EmptyDataError:
        # if no result in tmp_path
        return ''

    if not is_genome_fasta:
        # the fasta file is not whole genome fasta, get real genome position from it
        chrom_start = fimo_out['sequence_name'].apply(lambda i: i.split('-')[0])
        final_df = pd.DataFrame(chrom_start.str.split(':').tolist(),
                                columns=['sequence_name', 'region_start'])
        final_df['region_start'] = final_df['region_start'].astype(int)
        # Calculate the real genome position. fimo_out['start'] is relative to final_df['region_start']
        final_df['start'] = final_df['region_start'] + fimo_out['start']
        final_df['stop'] = final_df['region_start'] + fimo_out['stop']
        final_df['strand'] = fimo_out['strand']
        final_df['motif_id'] = fimo_out['motif_id']
        final_df['score'] = fimo_out['score']
        final_df['p-value'] = fimo_out['p-value']
    else:
        final_df = fimo_out

    # motif stats are not saved
    final_df = final_df.sort_values('score', ascending=False).iloc[:min(final_df.shape[0], top_n), :]
    print(f'{motif_path}, N motif total={final_df.shape[0]}')

    col_order = ['sequence_name', 'start', 'stop', 'motif_id', 'strand']
    final_df.loc[:, col_order].to_hdf(output_path, key='data')

    subprocess.run(['rm', '-f', tmp_path], check=True)
    return output_path


def _scan_motif_over_fasta(meme_motif_file, fasta_path, cpu, output_dir, path_to_fimo='',
                           raw_score_thresh=7, raw_p_value_thresh=5e-4, top_n=300000,
                           is_genome_fasta=False):
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
                                     is_genome_fasta=is_genome_fasta,
                                     top_n=top_n)
            futures[future] = (uid, name, motif_file_path)

        n = 0
        motif_records = []
        for future in as_completed(futures):
            output_path = future.result()
            uid, name, motif_file_path = futures[future]
            if output_path == '':
                # fimo return nothing, the parameter is too stringent for the motif.
                # by default parms, this rarely happens (2 in 1800).
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


def _aggregate_motif_beds(bed_file_paths, output_path, cpu=1, sort_mem_gbs=1):
    output_path = str(pathlib.Path(output_path).resolve())

    temp_bed = output_path + '.tmp_input.bed'
    with open(temp_bed, 'w') as temp_f:
        for bed_file_path in bed_file_paths:
            bed_df = pd.read_hdf(bed_file_path)[['sequence_name', 'start', 'stop', 'motif_id', 'strand']]
            csv_string = bed_df.to_csv(sep='\t', header=None, index=None)
            temp_f.write(csv_string)

    subprocess.run(shlex.split(f'sort -k1,1 -k2,2n --parallel={cpu} -S {sort_mem_gbs}G {temp_bed} -o {output_path}'),
                   stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', check=True)
    subprocess.run(shlex.split(f'rm -f {temp_bed}'))
    return output_path
