import pathlib
from collections import defaultdict

import msgpack
import pandas as pd

from .._bam_to_allc import _get_chromosome_sequence_upper, _read_faidx
from .._open import open_gz


def generate_c_motif_database(motif_bed, fasta_path):
    fai_path = f'{fasta_path}.fai'
    fai_df = _read_faidx(fai_path)

    cur_name_id = 0
    motif_names = {}
    with open_gz(motif_bed) as f:
        cur_chrom = ''
        full_chrom_seq = ''
        record_dict = defaultdict(list)

        for line in f:
            chrom, start, end, names, score, strand = line.strip('\n').split('\t')
            try:
                names = motif_names[names]
            except KeyError:
                motif_names[names] = cur_name_id
                names = cur_name_id
                cur_name_id += 1

            start = int(start)
            end = int(end)
            if chrom != cur_chrom:
                try:
                    full_chrom_seq = _get_chromosome_sequence_upper(fasta_path, fai_df, chrom)
                except KeyError:
                    continue
                # dump records and init
                if cur_chrom != '':
                    with open(f'{cur_chrom}.c_motif.msg', 'wb') as f:
                        msgpack.dump(record_dict, f)
                cur_chrom = chrom
                print(chrom)
                record_dict = defaultdict(list)

            motif = full_chrom_seq[start:end]
            if strand == '+':
                c_prefix = 1
                g_prefix = 0
            else:
                c_prefix = 0
                g_prefix = 1
            for i_base, base in enumerate(motif):
                pos = start + i_base + 1
                if base == 'C':
                    record_dict[pos].append((c_prefix, names))
                elif base == 'G':
                    record_dict[pos].append((g_prefix, names))

        # dump records and init
        if cur_chrom != '':
            with open(f'{cur_chrom}.c_motif.msg', 'wb') as f:
                msgpack.dump(record_dict, f)
    with open('motif_names.msg', 'wb') as f:
        msgpack.dump(motif_names, f)

    bin_size = 10000000

    # split msg into chrom bins
    for path in pathlib.Path().glob('chr*c_motif.msg'):
        bin_dict = defaultdict(dict)
        chrom = path.name.split('.')[0]
        print(chrom)
        with open(f'{chrom}.c_motif.msg', 'rb') as f:
            d = msgpack.unpackb(f.read(), raw=True, use_list=False)
            for k, v in d.items():
                nbin = k // bin_size
                bin_dict[nbin][k] = v
        for nbin, data in bin_dict.items():
            with open(f'human_hg19/c_motif/{chrom}.{int(nbin) * bin_size}-{int(nbin + 1) * bin_size}.c_motif.msg',
                      'wb') as f:
                msgpack.dump(data, f)

    # make lookup table
    c_motif_paths = list(pathlib.Path('./c').glob('chr*msg'))
    records = []
    for path in c_motif_paths:
        record = {'file_name': path.name,
                  'chrom': path.name.split('.')[0],
                  'start': path.name.split('.')[1].split('-')[0],
                  'end': path.name.split('.')[1].split('-')[1]}
        records.append(record)
    df = pd.DataFrame(records)[['chrom', 'start', 'end', 'file_name']]
    return
