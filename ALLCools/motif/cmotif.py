import pathlib
import subprocess
from collections import defaultdict

import msgpack
import pandas as pd

from .._doc import *
from .fimo import _scan_motif_over_fasta, _aggregate_motif_beds
from ..dmr.get_fasta import get_fasta
from .._bam_to_allc import _get_chromosome_sequence_upper, _read_faidx


def _generate_c_motif_database(motif_bed, fasta_path, output_dir, bin_size=10000000):
    fai_path = f'{fasta_path}.fai'
    fai_df = _read_faidx(fai_path)

    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir()

    cur_name_id = 0
    motif_names = {}
    chrom_file_paths = []
    with open(motif_bed) as f:
        cur_chrom = ''
        full_chrom_seq = ''
        record_dict = defaultdict(list)
        for line in f:
            chrom, start, end, name, strand = line.strip('\n').split('\t')
            try:
                motif_int_id = motif_names[name]
            except KeyError:
                motif_names[name] = cur_name_id
                motif_int_id = cur_name_id
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
                    file_path = output_dir / f'{cur_chrom}.c_motif.msg'
                    chrom_file_paths.append(file_path)
                    with open(file_path, 'wb') as chr_f:
                        msgpack.dump(record_dict, chr_f)
                cur_chrom = chrom
                print(chrom)
                record_dict = defaultdict(list)

            motif = full_chrom_seq[start:end]
            if strand == '+':
                c_prefix = 1  # for forward strand
                g_prefix = 0  # for reverse strand
            else:
                c_prefix = 0
                g_prefix = 1
            for i_base, base in enumerate(motif):
                pos = start + i_base + 1
                if base == 'C':
                    record_dict[pos].append((c_prefix, motif_int_id))
                elif base == 'G':
                    record_dict[pos].append((g_prefix, motif_int_id))

        # dump records and init
        if cur_chrom != '':
            file_path = output_dir / f'{cur_chrom}.c_motif.msg'
            chrom_file_paths.append(file_path)
            with open(file_path, 'wb') as chr_f:
                msgpack.dump(record_dict, chr_f)
    with open(output_dir / 'MOTIF_NAMES.msg', 'wb') as f:
        msgpack.dump(motif_names, f)

    # split msg into chrom bins
    bin_paths = []
    for path in chrom_file_paths:
        bin_dict = defaultdict(dict)
        chrom = path.name.split('.')[0]
        print(chrom)
        with open(path, 'rb') as f:
            d = msgpack.unpackb(f.read(), raw=True, use_list=False)
            for k, v in d.items():
                nbin = k // bin_size
                bin_dict[nbin][k] = v
        for nbin, data in bin_dict.items():
            bin_path = output_dir / f'{chrom}.{int(nbin) * bin_size}-{int(nbin + 1) * bin_size}.c_motif.msg'
            bin_paths.append(bin_path)
            with open(bin_path, 'wb') as f:
                msgpack.dump(data, f)
        subprocess.run(['rm', '-f', path])

    # make lookup table
    records = []
    for path in bin_paths:
        record = {'file_name': path.name,
                  'chrom': path.name.split('.')[0],
                  'start': int(path.name.split('.')[1].split('-')[0]),
                  'end': int(path.name.split('.')[1].split('-')[1])}
        records.append(record)
    df = pd.DataFrame(records)[['chrom', 'start', 'end', 'file_name']]
    df.to_hdf(output_dir / 'LOOKUP_TABLE.hdf', key='data')
    return


@doc_params(chrom_size_path_doc=chrom_size_path_doc,
            cpu_basic_doc=cpu_basic_doc)
def generate_cmotif_database(bed_file_paths,
                             reference_fasta,
                             motif_files,
                             output_dir,
                             slop_b=None,
                             chrom_size_path=None,
                             cpu=1,
                             sort_mem_gbs=5,
                             path_to_fimo='',
                             raw_score_thresh=8.,
                             raw_p_value_thresh=2e-4,
                             top_n=300000,
                             cmotif_bin_size=10000000):
    """\
    Generate lookup table for motifs all the cytosines belongs to.
    BED files are used to limit cytosine scan in certain regions.
    Scanning motif over whole genome is very noisy, better scan it in some functional part of genome.
    The result files will be in the output

    Parameters
    ----------
    bed_file_paths
        Paths of bed files. Multiple bed will be merged to get a final region set.
        The motif scan will only happen on the regions defined in these bed files.
    reference_fasta
        FASTA file path of the genome to scan
    motif_files
        MEME motif files that contains all the motif information.
    output_dir
        Output directory of C-Motif database
    slop_b
        Whether add slop to both ends of bed files.
    chrom_size_path
        {chrom_size_path_doc}
        Needed if slop_b is not None
    cpu
        {cpu_basic_doc}
    sort_mem_gbs
        Maximum memory usage in GBs when sort bed files
    path_to_fimo
        Path to fimo executable, if fimo is not in PATH
    raw_score_thresh
        Threshold of raw motif match likelihood score, see fimo doc for more info.
    raw_p_value_thresh
        Threshold of raw motif match P-value, see fimo doc for more info.
    top_n
        If too much motif found, will order them by likelihood score and keep top matches.
    cmotif_bin_size
        Bin size of single file in C-Motif database. No impact on results, better keep the default.
    Returns
    -------

    """
    try:
        if path_to_fimo != '':
            path_to_fimo = path_to_fimo.rstrip('/') + '/'
        subprocess.run([f'{path_to_fimo}fimo', '--version'], check=True)
    except subprocess.CalledProcessError as e:
        print('Test fimo failed. Make sure fimo (MEME suite) is installed and path_to_fimo set correctly.')
        raise e

    temp_paths = []

    output_dir = pathlib.Path(output_dir).absolute()
    bed_fasta_path = str(output_dir) + '.fa'
    temp_paths.append(bed_fasta_path)

    # step 1: Merge all bed files, get fasta file from the merged bed file.
    get_fasta(bed_file_paths=bed_file_paths,
              fasta_path=reference_fasta,
              output_path=bed_fasta_path,
              slop_b=slop_b,
              chrom_size_path=chrom_size_path,
              cpu=cpu, sort_mem_gbs=sort_mem_gbs)

    # step 2: scan motif over bed_fasta
    motif_scan_dir = str(output_dir) + '.motif_scan'
    temp_paths.append(motif_scan_dir)

    _scan_motif_over_fasta(
        meme_motif_file=motif_files,
        fasta_path=bed_fasta_path,
        cpu=cpu,
        output_dir=motif_scan_dir,
        path_to_fimo=path_to_fimo,
        raw_score_thresh=raw_score_thresh,
        raw_p_value_thresh=raw_p_value_thresh,
        top_n=top_n,
        is_genome_fasta=False)

    # step 3: aggregate motif bed
    total_motif_bed = str(output_dir) + '.total_motif.bed'
    lookup_table = pd.read_hdf(f'{motif_scan_dir}/LOOKUP_TABLE.hdf')
    bed_file_paths = lookup_table['output_file_path'].tolist()
    if len(bed_file_paths) == 0:
        print('None of the motif have any match in the current setting. Nothing to continue.'
              'Does the region to small or the threshold too stringent?')
        return
    _aggregate_motif_beds(bed_file_paths=bed_file_paths,
                          output_path=total_motif_bed,
                          cpu=cpu,
                          sort_mem_gbs=sort_mem_gbs)

    # step 4: generate cmotif dataset
    _generate_c_motif_database(motif_bed=total_motif_bed,
                               fasta_path=reference_fasta,
                               output_dir=output_dir,
                               bin_size=cmotif_bin_size)

    # step 5: cleaning
    for path in temp_paths:
        subprocess.run(['rm', '-rf', path])
    # not delete this file, but gzip it.
    subprocess.run(['bgzip', total_motif_bed])
    # save motif files into output_dir
    motif_file_dir = output_dir / 'MOTIF_FILE'
    motif_file_dir.mkdir(exist_ok=True)
    for motif_file in motif_files:
        new_path = motif_file_dir / pathlib.Path(motif_file).name
        subprocess.run(['cp', motif_file, new_path], check=True)
    return
