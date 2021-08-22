import collections
import functools
import gzip
import itertools
import os
import pathlib
import shlex
from collections import defaultdict
from functools import partial
from subprocess import run, PIPE, CalledProcessError
from typing import Union, List

import numpy as np
import pandas as pd

from ._doc import *
from ._open import open_allc

IUPAC_TABLE = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'R': 'AG',
    'Y': 'CT',
    'S': 'GC',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ATC',
    'V': 'ACG',
    'N': 'ATCGN'
}

COMPLIMENT_BASE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                   'a': 't', 'c': 'g', 't': 'a', 'g': 'c',
                   'N': 'N', 'n': 'n'}


def reverse_compliment(seq):
    return ''.join(map(lambda i: COMPLIMENT_BASE[i], seq[::-1]))


def get_allc_chroms(allc_path):
    p = run(['tabix', '-l', allc_path],
            check=True, stderr=PIPE, stdout=PIPE, encoding='utf8')
    return p.stdout.strip('\n').split('\n')


@functools.lru_cache(maxsize=100)
def parse_mc_pattern(pattern: str) -> set:
    """
    parse mC context pattern
    """
    # IUPAC DNA abbr. table
    all_pos_list = []
    pattern = pattern.upper()
    for base in pattern:
        try:
            all_pos_list.append(IUPAC_TABLE[base])
        except KeyError:
            raise KeyError(f'Base {base} is not in IUPAC table.')
    context_set = set([''.join(i) for i in itertools.product(*all_pos_list)])
    return context_set


@functools.lru_cache(maxsize=10)
def parse_chrom_size(path, remove_chr_list=None):
    """
    support simple UCSC chrom size file, or .fai format (1st and 2nd columns same as chrom size file)

    return chrom:length dict
    """
    if remove_chr_list is None:
        remove_chr_list = []

    with open(path) as f:
        chrom_dict = collections.OrderedDict()
        for line in f:
            # *_ for other format like fadix file
            chrom, length, *_ = line.strip('\n').split('\t')
            if chrom in remove_chr_list:
                continue
            chrom_dict[chrom] = int(length)
    return chrom_dict


def chrom_dict_to_id_index(chrom_dict, bin_size):
    sum_id = 0
    index_dict = {}
    for chrom, chrom_length in chrom_dict.items():
        index_dict[chrom] = sum_id
        sum_id += chrom_length // bin_size + 1
    return index_dict


def get_bin_id(chrom, chrom_index_dict, bin_start, bin_size) -> int:
    chrom_index_start = chrom_index_dict[chrom]
    n_bin = bin_start // bin_size
    return chrom_index_start + n_bin


def genome_region_chunks(chrom_size_path: str,
                         bin_length: int = 10000000,
                         combine_small: bool = True) -> List[str]:
    """
    Split the whole genome into bins, where each bin is {bin_length} bp. Used for tabix region query

    Parameters
    ----------
    chrom_size_path
        Path of UCSC genome size file
    bin_length
        length of each bin
    combine_small
        whether combine small regions into one record

    Returns
    -------
    list of records in tabix query format
    """
    chrom_size_dict = parse_chrom_size(chrom_size_path)

    cur_chrom_pos = 0
    records = []
    record_lengths = []
    for chrom, chrom_length in chrom_size_dict.items():
        while cur_chrom_pos + bin_length <= chrom_length:
            # tabix region is 1 based and inclusive
            records.append(f'{chrom}:{cur_chrom_pos}-{cur_chrom_pos + bin_length - 1}')
            cur_chrom_pos += bin_length
            record_lengths.append(bin_length)
        else:
            records.append(f'{chrom}:{cur_chrom_pos}-{chrom_length}')
            cur_chrom_pos = 0
            record_lengths.append(chrom_length - cur_chrom_pos)

    # merge small records (when bin larger then chrom length)
    final_records = []
    if combine_small:
        temp_records = []
        cum_length = 0
        for record, record_length in zip(records, record_lengths):
            temp_records.append(record)
            cum_length += record_length
            if cum_length >= bin_length:
                final_records.append(' '.join(temp_records))
                temp_records = []
                cum_length = 0
        if len(temp_records) != 0:
            final_records.append(' '.join(temp_records))
    else:
        for record in records:
            final_records.append(record)
    return final_records


def parse_file_paths(input_file_paths: Union[str, list]) -> list:
    if isinstance(input_file_paths, list) and (len(input_file_paths) == 1):
        input_file_paths = input_file_paths[0]

    if isinstance(input_file_paths, str):
        if '*' in input_file_paths:
            import glob
            file_list = glob.glob(input_file_paths)
        else:
            file_list = []
            with open(input_file_paths) as f:
                for line in f:
                    file_list.append(line.strip('\n'))
        _file_list = file_list
    elif isinstance(input_file_paths, list):
        _file_list = input_file_paths
    else:
        raise TypeError('File paths input is neither str nor list.')

    final_file_list = []
    for path in _file_list:
        real_path = pathlib.Path(path).resolve()
        if not real_path.exists():
            raise FileNotFoundError(f'{path} provided do not exist.')
        final_file_list.append(str(real_path))
    return _file_list


def get_md5(file_path):
    file_md5 = run(shlex.split(f'md5sum {file_path}'), stdout=PIPE, encoding='utf8', check=True).stdout
    file_md5 = file_md5.split(' ')[0]
    return file_md5


def check_tbi_chroms(file_path, genome_dict, same_order=False):
    file_tabix_path = pathlib.Path(str(file_path) + '.tbi')
    if not file_tabix_path.exists():
        print(f'{file_path} do not have .tbi index. Use tabix to index it.')
        return False

    tbi_time = os.path.getmtime(file_tabix_path)
    file_time = os.path.getmtime(file_path)
    if file_time > tbi_time:
        # tabix file create earlier than ALLC, something may changed for ALLC.
        return False

    try:
        chroms = run(['tabix', '-l', file_path],
                     stdout=PIPE,
                     encoding='utf8',
                     check=True).stdout.strip().split('\n')
        if len(set(chroms) - genome_dict.keys()) != 0:
            return False

        if same_order:
            ref_order = list(genome_dict.keys())
            cur_pos = -1
            for chrom in chroms:
                chrom_pos = ref_order.index(chrom)
                if chrom_pos < cur_pos:
                    return False
                else:
                    cur_pos = chrom_pos
    except CalledProcessError:
        return False
    return True


def generate_chrom_bin_bed_dataframe(chrom_size_path: str,
                                     window_size: int,
                                     step_size: int = None) -> pd.DataFrame:
    """
    Generate BED format dataframe based on UCSC chrom size file and window_size
    return dataframe contain 3 columns: chrom, start, end. The index is 0 based continue bin index.
    """
    if step_size is None:
        step_size = window_size
    chrom_size_dict = parse_chrom_size(chrom_size_path)
    records = []
    for chrom, chrom_length in chrom_size_dict.items():
        bin_start = np.array(list(range(0, chrom_length, step_size)))
        bin_end = bin_start + window_size
        bin_end[np.where(bin_end > chrom_length)] = chrom_length
        chrom_df = pd.DataFrame(dict(bin_start=bin_start, bin_end=bin_end))
        chrom_df['chrom'] = chrom
        records.append(chrom_df)
    total_df = pd.concat(records)[['chrom', 'bin_start', 'bin_end']].reset_index(drop=True)
    return total_df


@doc_params(allc_path_doc=allc_path_doc)
def profile_allc(allc_path, drop_n=True, n_rows=1000000, output_path=None):
    """\
    Generate some summary statistics of 1 ALLC.
    1e8 rows finish in about 5 min.

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    drop_n
        Whether to drop context that contain N, such as CCN.
        This is usually very rare and need to be dropped.
    n_rows
        Number of rows to calculate the profile from.
        The default number is usually sufficient to get pretty precise assumption.
    output_path
        Path of the output file. If None, will save the profile next to input ALLC file.
    Returns
    -------

    """
    # TODO write test
    # Find best default value
    if 'gz' in allc_path:
        opener = partial(gzip.open, mode='rt')
    else:
        opener = partial(open, mode='r')

    # initialize count dict
    mc_sum_dict = defaultdict(int)
    cov_sum_dict = defaultdict(int)
    cov_sum2_dict = defaultdict(int)  # sum of square, for calculating variance
    rate_sum_dict = defaultdict(float)
    rate_sum2_dict = defaultdict(float)  # sum of square, for calculating variance
    context_count_dict = defaultdict(int)
    with opener(allc_path) as f:
        n = 0
        for line in f:
            chrom, pos, strand, context, mc, cov, p = line.split('\t')
            if drop_n and 'N' in context:
                continue
            # mc and cov
            mc_sum_dict[context] += int(mc)
            cov_sum_dict[context] += int(cov)
            cov_sum2_dict[context] += int(cov) ** 2
            # raw base rate
            rate = int(mc) / int(cov)
            rate_sum_dict[context] += rate
            rate_sum2_dict[context] += rate ** 2
            # count context finally
            context_count_dict[context] += 1
            n += 1
            if (n_rows is not None) and (n >= n_rows):
                break

    # overall count
    profile_df = pd.DataFrame({'partial_mc': mc_sum_dict,
                               'partial_cov': cov_sum_dict})
    profile_df['base_count'] = pd.Series(context_count_dict)
    profile_df['overall_mc_rate'] = profile_df['partial_mc'] / profile_df['partial_cov']

    # cov base mean and base std.
    # assume that base cov follows normal distribution
    cov_sum_series = pd.Series(cov_sum_dict)
    cov_sum2_series = pd.Series(cov_sum2_dict)
    profile_df['base_cov_mean'] = cov_sum_series / profile_df['base_count']
    profile_df['base_cov_std'] = np.sqrt(
        (cov_sum2_series / profile_df['base_count']) - profile_df['base_cov_mean'] ** 2)

    # assume that base rate follow beta distribution
    # so that observed rate actually follow joint distribution of beta (rate) and normal (cov) distribution
    # here we use the observed base_rate_mean and base_rate_var to calculate
    # approximate alpha and beta value for the base rate beta distribution
    rate_sum_series = pd.Series(rate_sum_dict)
    rate_sum2_series = pd.Series(rate_sum2_dict)
    profile_df['base_rate_mean'] = rate_sum_series / profile_df['base_count']
    profile_df['base_rate_var'] = (rate_sum2_series / profile_df['base_count']) - profile_df['base_rate_mean'] ** 2

    # based on beta distribution mean, var
    # a / (a + b) = base_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = base_rate_var
    # we have:
    a = (1 - profile_df['base_rate_mean']) * (profile_df['base_rate_mean'] ** 2) / profile_df['base_rate_var'] - \
        profile_df['base_rate_mean']
    b = a * (1 / profile_df['base_rate_mean'] - 1)
    profile_df['base_beta_a'] = a
    profile_df['base_beta_b'] = b

    if output_path is not None:
        profile_df.to_csv(output_path, sep='\t')
        return None
    else:
        return profile_df


@doc_params(allc_path_doc=allc_path_doc)
def tabix_allc(allc_path, reindex=False):
    """\
    a simple wrapper of tabix command to index 1 ALLC file

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    reindex
        If True, will force regenerate the ALLC index.
    Returns
    -------

    """
    if os.path.exists(f'{allc_path}.tbi') and not reindex:
        return
    run(shlex.split(f'tabix -f -b 2 -e 2 -s 1 {allc_path}'),
        check=True)
    return


@doc_params(allc_path_doc=allc_path_doc,
            chrom_size_path_doc=chrom_size_path_doc,
            compress_level_doc=compress_level_doc,
            remove_additional_chrom_doc=remove_additional_chrom_doc)
def standardize_allc(allc_path, chrom_size_path, compress_level=5,
                     remove_additional_chrom=False):
    """\
    Standardize 1 ALLC file by checking:
        1. No header in the ALLC file;
        2. Chromosome names in ALLC must be same as those in the chrom_size_path file, including "chr";
        3. Output file will be bgzipped with .tbi index
        4. Remove additional chromosome (remove_additional_chrom=True) or
           raise KeyError if unknown chromosome found (default)

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    chrom_size_path
        {chrom_size_path_doc}
    compress_level
        {compress_level_doc}
    remove_additional_chrom
        {remove_additional_chrom_doc}
    Returns
    -------

    """
    genome_dict = parse_chrom_size(chrom_size_path)
    if check_tbi_chroms(allc_path, genome_dict):
        # means ALLC is already standard
        return

    if 'chr1' in genome_dict:
        raw_add_chr = True
    else:
        raw_add_chr = False
    with open_allc(allc_path) as f, \
            open_allc(allc_path + '.tmp.gz', mode='w',
                      compresslevel=compress_level) as wf:
        cur_chrom = "TOTALLY_NOT_A_CHROM"
        cur_start = cur_chrom + '\t'
        cur_pointer = 0
        index_lines = []
        buffer_lines = ''
        line_count = 0
        add_chr = raw_add_chr
        for line in f:
            if line_count == 0:
                # for very old allc files, which contain header line
                ll = line.split('\t')
                try:
                    int(ll[1])  # pos
                    int(ll[4])  # mc
                    int(ll[5])  # cov
                except ValueError:
                    # The first line is header, remove header
                    continue
            if line_count < 2:
                # 1st line could be header that startswith chr
                # so check 1st and 2nd row
                if line.startswith('chr'):
                    add_chr = False
            if add_chr:
                line = 'chr' + line
            if not line.startswith(cur_start):
                fields = line.split("\t")
                cur_chrom = fields[0]
                if cur_chrom not in genome_dict:
                    if not remove_additional_chrom:
                        raise KeyError(f'{cur_chrom} not exist in genome size file, '
                                       f'set remove_additional_chrom=True if want to remove additional chroms')
                    else:
                        # skip unknown chromosome
                        continue
                index_lines.append(cur_chrom + "\t" + str(cur_pointer) + "\n")
                cur_start = cur_chrom + '\t'
            cur_pointer += len(line)
            buffer_lines += line
            line_count += 1
            if line_count % 50000 == 0:
                wf.write(buffer_lines)
                buffer_lines = ''
        wf.write(buffer_lines)

    run(shlex.split(f'mv {allc_path} {allc_path}.bp'), check=True)
    run(shlex.split(f'mv {allc_path}.tmp.gz {allc_path}'), check=True)
    run(shlex.split(f'rm -f {allc_path}.bp'), check=True)
    tabix_allc(allc_path, reindex=True)
    return


def _transfer_bin_size(bin_size: int) -> str:
    """Get proper str for a large bin_size"""
    if bin_size > 1000000:
        bin_size_mode = bin_size % 1000000
        bin_size_mode = f'{bin_size_mode / 1000000:.1f}'[1:] if bin_size_mode >= 100000 else ''
        bin_size_str = f'{bin_size // 1000000}{bin_size_mode}m'
    elif bin_size > 1000:
        bin_size_mode = bin_size % 1000
        bin_size_mode = f'{bin_size_mode / 1000:.1f}'[1:] if bin_size_mode >= 100 else ''
        bin_size_str = f'{bin_size // 1000}{bin_size_mode}k'
    else:
        bin_size_str = f'{bin_size}'
    return bin_size_str


def parse_dtype(dtype):
    if isinstance(dtype, str):
        if dtype == 'uint8':
            return np.uint16
        elif dtype == 'uint16':
            return np.uint16
        elif dtype == 'uint32':
            return np.uint32
        elif dtype == 'uint64':
            return np.uint64
        elif dtype == 'int8':
            return np.int8
        elif dtype == 'int16':
            return np.int16
        elif dtype == 'int32':
            return np.int32
        elif dtype == 'int64':
            return np.int64
        elif dtype == 'bool':
            return np.bool
        else:
            raise ValueError(f'Unknown dtype {dtype}')
    else:
        return dtype


def binary_count(mc, cov):
    if mc == 0:
        return 0, 1
    elif mc == cov:
        return 1, 1
    else:
        return 0, 0
