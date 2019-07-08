"""
Run FIMO to scan motif over FASTA sequences
"""
import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd


def _split_meme_motif_file(meme_motif_path, output_dir):
    """
    Given multi motif meme format file, split into single motif meme format file

    Parameters
    ----------
    meme_motif_path
    output_dir

    Returns
    -------

    """
    with open(meme_motif_path) as f:
        first_line = f.readline()
        if not first_line.startswith('MEME version'):
            raise ValueError('Input file need to be MEME motif format.')

        f.seek(0)
        header = True
        records = {}
        header_text = ''
        motif_tmp_text = ''
        for line in f:
            if line.startswith('MOTIF'):
                try:
                    _, uid, name = line.strip('\n').split(' ')
                except ValueError:
                    _, uid = line.strip('\n').split(' ')
                    name = ''
                header = False

                if motif_tmp_text != '':
                    records[(uid, name)] = header_text + motif_tmp_text
                motif_tmp_text = line
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


def _fimo_runner(motif_path, fasta_path, output_path):
    """Run fimo for a single motif over single fasta file"""
    cmd = f'fimo --text --skip-matched-sequence {motif_path} {fasta_path}'

    tmp_path = str(output_path) + "tmp"
    with open(tmp_path, 'w') as out:
        try:
            subprocess.run(shlex.split(cmd),
                           stdout=out,
                           stderr=subprocess.PIPE,
                           encoding='utf8',
                           check=True)
        except subprocess.CalledProcessError as e:
            print(cmd)
            print(e.stdout)
            raise e

    col_order = ['sequence_name', 'start', 'stop',
                 'motif_id', 'motif_alt_id', 'strand',
                 'score', 'p-value']
    try:
        fimo_out = pd.read_csv(tmp_path, sep='\t')
        # same as sort -k1,1 -k2,2n, so in bedtools intersect, we can use --sorted
        fimo_out = fimo_out.sort_values(by=['sequence_name', 'start']).reset_index(drop=True)
    except pd.errors.EmptyDataError:
        # if no result in tmp_path
        fimo_out = pd.DataFrame([], columns=col_order)

    fimo_out.loc[:, col_order].to_msgpack(output_path, compress='zlib')
    subprocess.run(['rm', '-f', tmp_path], check=True)
    return output_path


def scan_motif_over_fasta(meme_motif_file, fasta_path, cpu, output_dir):
    output_dir = pathlib.Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    motif_file_records = _split_meme_motif_file(meme_motif_file, output_dir)
    print(f'{len(motif_file_records)} motifs to count.')

    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for uid, name, motif_file_path in motif_file_records:
            future = executor.submit(_fimo_runner,
                                     motif_path=motif_file_path,
                                     fasta_path=fasta_path,
                                     output_path=output_dir / (pathlib.Path(motif_file_path).name[:-5] + '.bed.msg'))
            futures[future] = (uid, name, motif_file_path)

        n = 0
        motif_records = []
        for future in as_completed(futures):
            output_path = future.result()
            uid, name, motif_file_path = futures[future]
            motif_records.append((uid, name, motif_file_path, output_path))
            n += 1
            if n % 100 == 0:
                print(f'Finish count {n} motifs.')

    look_up_table = pd.DataFrame(motif_records,
                                 columns=['uid', 'name', 'motif_file_path', 'output_file_path']).set_index('uid')
    look_up_table.to_csv(output_dir / 'LOOK_UP_TABLE.tsv', sep='\t')
    return
