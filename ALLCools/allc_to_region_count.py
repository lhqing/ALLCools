import pathlib
import shlex
import subprocess
from collections import defaultdict
from functools import partial
from typing import Union, Tuple

from ._open_ import open_allc, open_gz
from .extract_allc import extract_allc
from .utilities import check_tbi_chroms, parse_chrom_size, parse_mc_pattern


def _allc_to_site_bed(allc_path: str,
                      out_prefix: str,
                      chrom_size_path: str,
                      mc_contexts: Union[str, list],
                      max_cov_cutoff: int = 9999) -> Tuple[list, list]:
    """
    1. Split allc by context into several site-bed files
    2. Keep chromosome sorted as the chrom_size_path
    3. BED file 4th column is mc, 5th column is cov, and each row is a single site
    """
    # TODO add parallel to this
    out_prefix = out_prefix.rstrip('.')
    if isinstance(mc_contexts, str):
        mc_contexts = mc_contexts.split(' ')
    mc_contexts = list(set(mc_contexts))

    # because mc_contexts can overlap (e.g. CHN, CAN)
    # each context may associate to multiple handle
    context_handle = defaultdict(list)
    handle_collect = []
    out_paths = []
    for mc_context in mc_contexts:
        out_path = out_prefix + f'.extract_{mc_context}.bed.gz'
        out_paths.append(out_path)
        _handle = open_allc(out_path, 'w')
        handle_collect.append(_handle)
        parsed_context_set = parse_mc_pattern(mc_context)
        for pattern in parsed_context_set:
            context_handle[pattern].append(_handle)

    # split file first
    chrom_size_dict = parse_chrom_size(chrom_size_path)
    with open_allc(allc_path, region=' '.join(chrom_size_dict.keys())) as allc:
        for line in allc:
            chrom, pos, _, context, mc, cov, *_ = line.strip('\n').split('\t')
            if int(cov) > max_cov_cutoff:
                continue
            bed_line = '\t'.join([chrom, pos, pos, mc, cov]) + '\n'
            try:
                [h.write(bed_line) for h in context_handle[context]]
            except KeyError:
                continue
    for handle in handle_collect:
        handle.close()
    return mc_contexts, out_paths


def _bedtools_map(region_bed, site_bed, out_bed, save_zero_cov=True):
    cmd = f'bedtools map -a {region_bed} -b {site_bed} -c 4,5 -o sum,sum'
    bed_out = subprocess.run(shlex.split(cmd),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             encoding='utf8', check=True)
    if out_bed.endswith('gz'):
        opener = partial(open_gz, mode='wt')
    else:
        opener = partial(open, mode='w')

    with opener(out_bed) as out_handle:
        for line in bed_out.stdout:
            if save_zero_cov:
                out_handle.write(line)
            else:
                # the last item is cov
                if not line.endswith('\t.\n'):
                    out_handle.write(line)
    return


def _chrom_dict_to_id_index(chrom_dict, bin_size):
    sum_id = 0
    index_dict = {}
    for chrom, chrom_length in chrom_dict.items():
        index_dict[chrom] = sum_id
        sum_id += chrom_length // bin_size + 1
        if chrom_length % bin_size == 0:
            sum_id += 1
    return index_dict


def _get_bin_id(chrom, chrom_index_dict, bin_start, bin_size) -> int:
    chrom_index_start = chrom_index_dict[chrom]
    n_bin = bin_start // bin_size
    return chrom_index_start + n_bin


def _map_to_sparse_chrom_bin(input_path, output_path, chrom_size_file,
                             bin_size=500):
    """
    Calculate chromosome bins regional count, output is sparse,
    bin_id constructed from chrom_size_file and can be reproduce.
    """
    chrom_dict = parse_chrom_size(chrom_size_file)
    chrom_index_dict = _chrom_dict_to_id_index(chrom_dict, bin_size)
    cur_chrom = 'TOTALLY_NOT_A_CHROM'
    cur_chrom_end = 0
    bin_end = min(cur_chrom_end, bin_size)
    bin_start = 0
    bin_id = -1
    temp_mc, temp_cov = -1, -1

    with open_gz(input_path) as in_handle, \
            open_gz(output_path, 'w') as out_handle:
        # add header to indicate chromosome order
        for line in in_handle:
            # site-bed format
            chrom, pos, _, mc, cov = line.split("\t")
            pos = int(pos)
            mc = int(mc)
            cov = int(cov)
            if pos >= bin_end or cur_chrom != chrom:
                # write line
                if temp_cov > 0:
                    out_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end,
                                                         bin_id, temp_mc, temp_cov])) + "\n")

                # reset_chrom
                if cur_chrom != chrom:
                    cur_chrom = chrom
                    cur_chrom_end = chrom_dict[chrom]

                # reset bin
                temp_mc, temp_cov = mc, cov
                bin_start = int(pos // bin_size * bin_size)
                bin_end = min(cur_chrom_end, bin_start + bin_size)
                bin_id = _get_bin_id(cur_chrom, chrom_index_dict, bin_start, bin_size)
            else:
                temp_mc += mc
                temp_cov += cov
        # write last piece
        if temp_cov > 0:
            out_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end,
                                                 bin_id, temp_mc, temp_cov])) + "\n")
    return output_path


def _transfer_bin_size(bin_size: int) -> str:
    """Get proper str for a large bin_size"""
    if bin_size > 1000000:
        bin_size_mode = bin_size % 1000000
        bin_size_mode = f'{bin_size_mode/1000000:.1f}'[1:] if bin_size_mode >= 100000 else ''
        bin_size_str = f'{bin_size//1000000}{bin_size_mode}m'
    elif bin_size > 1000:
        bin_size_mode = bin_size % 1000
        bin_size_mode = f'{bin_size_mode/1000:.1f}'[1:] if bin_size_mode >= 100 else ''
        bin_size_str = f'{bin_size//1000}{bin_size_mode}k'
    else:
        bin_size_str = f'{bin_size}'
    return bin_size_str


def allc_to_region_count(allc_path,
                         output_prefix,
                         chrom_size_path,
                         mc_contexts,
                         split_strand=True,
                         region_bed_paths=None,
                         region_bed_names=None,
                         bin_sizes=None,
                         max_cov_cutoff=9999,
                         save_zero_cov=True,
                         remove_tmp=True):
    """\
    Calculate mC and cov at regional level. Region accepted in 2 forms:
    1. BED file, provided by region_bed_paths, containing arbitrary regions and use bedtools map to calculate
    2. Fix-size non-overlap genome bins, provided by bin_sizes, this is much faster to calculate than 1.
    The output is in 6-column bed-like format:
    chrom   start   end region_uid  mc  cov

    Parameters
    ----------
    allc_path
        Path to the ALLC file
    output_prefix
        Path to output prefix
    chrom_size_path
        Path to UCSC chrom size file
    mc_contexts
        mC context list to calculate
    region_bed_paths
        Path to BED files
    region_bed_names
        Matched name for each BED file provided in region_bed_paths
    bin_sizes
        Sizes of genome bins to calculate
    max_cov_cutoff
        Max cov filter for a single site in ALLC
    save_zero_cov
        Whether to save the regions that have 0 cov, only apply to region count but not the chromosome count
    remove_tmp
        Whether to remove the temporary file
    Returns
    -------
    """
    genome_dict = parse_chrom_size(chrom_size_path)
    if bin_sizes is None and region_bed_paths is None:
        raise ValueError('Either bin_sizes or region_bed_paths should be provided.')

    # check bed file
    # 1. bgzip and tabix
    # 2. order of chrom should be the same as order of chrom_size_path
    if len(region_bed_names) != len(region_bed_paths):
        raise ValueError('Different number of bed names and paths')
    for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
        if not check_tbi_chroms(region_bed_path, genome_dict):
            raise ValueError(f'The bed file {region_bed_path} chromosome order is different '
                             f'from the {chrom_size_path}')

    print('Extract ALLC context')
    output_prefix = output_prefix.rstrip('.')
    strandness = 'split' if split_strand else 'both'
    output_paths = extract_allc(allc_path=allc_path,
                                output_prefix=output_prefix,
                                mc_contexts=mc_contexts,
                                strandness=strandness,
                                output_format='bed5',
                                region=None,
                                cov_cutoff=max_cov_cutoff)
    path_dict = {}
    for path in output_paths:
        # this is according to extract_allc name pattern
        info_type = pathlib.Path(path).name.split('.')[-4]  # {mc_context}-{strandness}
        path_dict[info_type] = path

    if region_bed_paths is not None:
        print('Map to regions.')
    save_flag = 'full' if save_zero_cov else 'sparse'
    for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
        for info_type, site_bed_path in path_dict.items():
            try:
                _bedtools_map(region_bed=region_bed_path,
                              site_bed=site_bed_path,
                              out_bed=output_prefix + f'.{region_name}_{info_type}.{save_flag}.bed.gz',
                              save_zero_cov=save_zero_cov)
            except subprocess.CalledProcessError as e:
                print(e.stderr)
                raise e

    if bin_sizes is not None:
        print('Map to chromosome bins.')
    for bin_size in bin_sizes:
        for info_type, site_bed_path in path_dict.items():
            _map_to_sparse_chrom_bin(input_path=site_bed_path,
                                     output_path=output_prefix + f'.chrom{_transfer_bin_size(bin_size)}'
                                                                 f'_{info_type}.sparse.bed.gz',
                                     chrom_size_file=chrom_size_path,
                                     bin_size=bin_size)

    if remove_tmp:
        print('Remove temporal files.')
        for site_bed_path in path_dict.values():
            subprocess.run(['rm', '-f', site_bed_path])
    return


def aggregate_region_count_matrix(region_count_table, region_count_tables, region_name, region_bed_path, bin_size):
    # TODO
    # this should only deal with a simple case, aggregate 2 sample*feature 2-D matrix, one for mc, one for cov,
    # output to full or sparse format
    return
