"""
This file is modified from methylpy https://github.com/yupenghe/methylpy

Author: Yupeng He
"""

import logging
import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import collections
import pandas as pd

from ._doc import *
from ._open import open_allc, open_bam
from .utilities import genome_region_chunks

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def _read_faidx(faidx_path):
    """
    Read fadix of reference fasta file
    samtools fadix ref.fa
    """
    return pd.read_csv(faidx_path, index_col=0, header=None, sep='\t',
                       names=['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH'])


def _get_chromosome_sequence_upper(fasta_path, fai_df, query_chrom):
    """
    read a whole chromosome sequence into memory
    """
    chrom_pointer = fai_df.loc[query_chrom, 'OFFSET']
    tail = fai_df.loc[query_chrom, 'LINEBASES'] - fai_df.loc[query_chrom, 'LINEWIDTH']
    seq = ""
    with open(fasta_path, 'r') as f:
        f.seek(chrom_pointer)
        for line in f:
            if line[0] == '>':
                break
            seq += line[:tail]  # trim \n
    return seq.upper()


def _get_bam_chrom_index(bam_path):
    result = subprocess.run(['samtools', 'idxstats', bam_path],
                            stdout=subprocess.PIPE, encoding='utf8').stdout

    chrom_set = set()
    for line in result.split('\n'):
        chrom = line.split('\t')[0]
        if chrom not in ['', '*']:
            chrom_set.add(chrom)
    return pd.Index(chrom_set)


def _bam_to_allc_worker(bam_path, reference_fasta, fai_df, output_path, region=None,
                        num_upstr_bases=0, num_downstr_bases=2,
                        buffer_line_number=100000, min_mapq=0, min_base_quality=1,
                        compress_level=5, tabix=True, save_count_df=False):
    """
    None parallel bam_to_allc worker function, call by bam_to_allc
    """
    # mpileup
    if region is None:
        mpileup_cmd = f"samtools mpileup -Q {min_base_quality} " \
                      f"-q {min_mapq} -B -f {reference_fasta} {bam_path}"
        pipes = subprocess.Popen(shlex.split(mpileup_cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
    else:
        bam_handle = open_bam(bam_path, region=region, mode='r',
                              include_header=True, samtools_parms_str=None)
        mpileup_cmd = f"samtools mpileup -Q {min_base_quality} " \
                      f"-q {min_mapq} -B -f {reference_fasta} -"
        pipes = subprocess.Popen(shlex.split(mpileup_cmd),
                                 stdin=bam_handle.file,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)

    result_handle = pipes.stdout

    output_file_handler = open_allc(output_path, mode='w',
                                    compresslevel=compress_level)

    # initialize variables
    complement = {"A": "T",
                  "T": "A",
                  "C": "G",
                  "G": "C",
                  "N": "N"}
    mc_sites = {'C', 'G'}
    context_len = num_upstr_bases + 1 + num_downstr_bases
    cur_chrom = ""
    line_counts = 0
    total_line = 0
    out = ""
    seq = None  # whole cur_chrom seq
    chr_out_pos_list = []
    cur_out_pos = 0
    cov_dict = collections.defaultdict(int)  # context: cov_total
    mc_dict = collections.defaultdict(int)  # context: mc_total

    # process mpileup result
    for line in result_handle:
        total_line += 1
        fields = line.split("\t")
        fields[2] = fields[2].upper()
        # if chrom changed, read whole chrom seq from fasta
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            chr_out_pos_list.append((cur_chrom, str(cur_out_pos)))
            # get seq for cur_chrom
            seq = _get_chromosome_sequence_upper(reference_fasta, fai_df, cur_chrom)

        if fields[2] not in mc_sites:
            continue

        # indels
        read_bases = fields[4]
        incons_basecalls = read_bases.count("+") + read_bases.count("-")
        if incons_basecalls > 0:
            read_bases_no_indel = ""
            index = 0
            prev_index = 0
            while index < len(read_bases):
                if read_bases[index] == "+" or read_bases[index] == "-":
                    # get insert size
                    indel_size = ""
                    ind = index + 1
                    while True:
                        try:
                            int(read_bases[ind])
                            indel_size += read_bases[ind]
                            ind += 1
                        except:
                            break
                    try:
                        # sometimes +/- does not follow by a number and
                        # it should be ignored
                        indel_size = int(indel_size)
                    except:
                        index += 1
                        continue
                    read_bases_no_indel += read_bases[prev_index:index]
                    index = ind + indel_size
                    prev_index = index
                else:
                    index += 1
            read_bases_no_indel += read_bases[prev_index:index]
            fields[4] = read_bases_no_indel

        # count converted and unconverted bases
        if fields[2] == "C":
            pos = int(fields[1]) - 1
            try:
                context = seq[(pos - num_upstr_bases):(pos + num_downstr_bases + 1)]
            except:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(".")
            converted_c = fields[4].count("T")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = "\t".join([cur_chrom, str(pos + 1), "+", context,
                                  str(unconverted_c), str(cov), "1"]) + "\n"
                cov_dict[context] += cov
                mc_dict[context] += unconverted_c
                out += data
                cur_out_pos += len(data)

        elif fields[2] == "G":
            pos = int(fields[1]) - 1
            try:
                context = "".join([complement[base]
                                   for base in reversed(
                        seq[(pos - num_downstr_bases):(pos + num_upstr_bases + 1)]
                    )]
                                  )
            except:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(",")
            converted_c = fields[4].count("a")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = "\t".join([cur_chrom, str(pos + 1), "-", context,
                                  str(unconverted_c), str(cov), "1"]) + "\n"
                cov_dict[context] += cov
                mc_dict[context] += unconverted_c
                out += data
                cur_out_pos += len(data)

        if line_counts > buffer_line_number:
            output_file_handler.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_file_handler.write(out)
    result_handle.close()
    output_file_handler.close()

    if tabix:
        subprocess.run(shlex.split(f'tabix -b 2 -e 2 -s 1 {output_path}'),
                       check=True)

    count_df = pd.DataFrame({'mc': mc_dict, 'cov': cov_dict})
    count_df['mc_rate'] = count_df['mc'] / count_df['cov']

    total_genome_length = fai_df['LENGTH'].sum()
    count_df['genome_cov'] = total_line / total_genome_length

    if save_count_df:
        count_df.to_csv(output_path + '.count.csv')
        return None
    else:
        return count_df


def _aggregate_count_df(count_dfs):
    total_df = pd.concat(count_dfs)
    total_df = total_df.groupby(total_df.index).sum()
    total_df['mc_rate'] = total_df['mc'] / total_df['cov']
    total_df['mc'] = total_df['mc'].astype(int)
    total_df['cov'] = total_df['cov'].astype(int)
    return total_df


@doc_params(compress_level_doc=compress_level_doc,
            cpu_basic_doc=cpu_basic_doc,
            reference_fasta_doc=reference_fasta_doc)
def bam_to_allc(bam_path,
                reference_fasta,
                output_path=None,
                cpu=1,
                num_upstr_bases=0,
                num_downstr_bases=2,
                min_mapq=10,
                min_base_quality=20,
                compress_level=5,
                save_count_df=False):
    """\
    Generate 1 ALLC file from 1 position sorted BAM file via samtools mpileup.

    Parameters
    ----------
    bam_path
        Path to 1 position sorted BAM file
    reference_fasta
        {reference_fasta_doc}
    output_path
        Path to 1 output ALLC file
    cpu
        {cpu_basic_doc} DO NOT use cpu > 1 for single cell ALLC generation.
        Parallel on cell level is better for single cell project.
    num_upstr_bases
        Number of upstream base(s) of the C base to include in ALLC context column,
        usually use 0 for normal BS-seq, 1 for NOMe-seq.
    num_downstr_bases
        Number of downstream base(s) of the C base to include in ALLC context column,
        usually use 2 for both BS-seq and NOMe-seq.
    min_mapq
        Minimum MAPQ for a read being considered, samtools mpileup parameter, see samtools documentation.
    min_base_quality
        Minimum base quality for a base being considered, samtools mpileup parameter,
        see samtools documentation.
    compress_level
        {compress_level_doc}
    save_count_df
        If true, save an ALLC context count table next to ALLC file.

    Returns
    -------
    count_df
        a pandas.DataFrame for overall mC and cov count separated by mC context.
    """
    buffer_line_number = 100000
    tabix = True

    # Check fasta index
    if not pathlib.Path(reference_fasta + ".fai").exists():
        raise FileNotFoundError("Reference fasta not indexed. Use samtools faidx to index it and run again.")
    fai_df = _read_faidx(pathlib.Path(reference_fasta + ".fai"))

    if not pathlib.Path(bam_path + ".bai").exists():
        subprocess.check_call(shlex.split("samtools index " + bam_path))

    # check chromosome between BAM and FASTA
    # samtools have a bug when chromosome not match...
    bam_chroms_index = _get_bam_chrom_index(bam_path)
    unknown_chroms = [i for i in bam_chroms_index if i not in fai_df.index]
    if len(unknown_chroms) != 0:
        unknown_chroms = ' '.join(unknown_chroms)
        raise IndexError(f'BAM file contain unknown chromosomes: {unknown_chroms}\n'
                         'Make sure you use the same genome FASTA file for mapping and bam-to-allc.')

    # if parallel, chunk genome
    if cpu > 1:
        regions = genome_region_chunks(reference_fasta + ".fai",
                                       bin_length=100000000,
                                       combine_small=False)
    else:
        regions = None

    # Output path
    input_path = pathlib.Path(bam_path)
    file_dir = input_path.parent
    if output_path is None:
        allc_name = 'allc_' + input_path.name.split('.')[0] + '.tsv.gz'
        output_path = str(file_dir / allc_name)
    else:
        if not output_path.endswith('.gz'):
            output_path += '.gz'

    if cpu > 1:
        raise NotImplementedError

        temp_out_paths = []
        for batch_id, region in enumerate(regions):
            temp_out_paths.append(output_path + f'.batch_{batch_id}.tmp.tsv.gz')

        with ProcessPoolExecutor(max_workers=cpu) as executor:
            future_dict = {}
            for batch_id, (region, out_temp_path) in enumerate(zip(regions, temp_out_paths)):
                _kwargs = dict(bam_path=bam_path,
                               reference_fasta=reference_fasta,
                               fai_df=fai_df,
                               output_path=out_temp_path,
                               region=region,
                               num_upstr_bases=num_upstr_bases,
                               num_downstr_bases=num_downstr_bases,
                               buffer_line_number=buffer_line_number,
                               min_mapq=min_mapq,
                               min_base_quality=min_base_quality,
                               compress_level=compress_level,
                               tabix=False,
                               save_count_df=False)
                future_dict[executor.submit(_bam_to_allc_worker, **_kwargs)] = batch_id

            count_dfs = []
            for future in as_completed(future_dict):
                batch_id = future_dict[future]
                try:
                    count_dfs.append(future.result())
                except Exception as exc:
                    log.info('%r generated an exception: %s' % (batch_id, exc))

            # aggregate ALLC
            with open_allc(output_path, mode='w', compresslevel=compress_level,
                           threads=1, region=None) as out_f:
                # TODO: Parallel ALLC is overlaped, the split by region strategy only split reads, but reads can overlap
                # need to adjust and merge end rows in aggregate ALLC
                for out_temp_path in temp_out_paths:
                    with open_allc(out_temp_path) as f:
                        for line in f:
                            out_f.write(line)
                    subprocess.check_call(['rm', '-f', out_temp_path])

            # tabix
            if tabix:
                subprocess.run(shlex.split(f'tabix -b 2 -e 2 -s 1 {output_path}'),
                               check=True)

            # aggregate count_df
            count_df = _aggregate_count_df(count_dfs)
            if save_count_df:
                count_df.to_csv(output_path + '.count.csv')
            return count_df
    else:
        return _bam_to_allc_worker(bam_path, reference_fasta, fai_df, output_path,
                                   region=None,
                                   num_upstr_bases=num_upstr_bases,
                                   num_downstr_bases=num_downstr_bases,
                                   buffer_line_number=buffer_line_number,
                                   min_mapq=min_mapq,
                                   min_base_quality=min_base_quality,
                                   compress_level=compress_level, tabix=tabix,
                                   save_count_df=save_count_df)
