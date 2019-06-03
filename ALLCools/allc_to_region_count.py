import pathlib
import shlex
import subprocess
from functools import partial
from typing import List

from ._doc_ import *
from ._open_ import open_gz
from .extract_allc import extract_allc
from .utilities import check_tbi_chroms, parse_chrom_size


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


@doc_params(allc_path_doc=allc_path_doc, chrom_size_path_doc=chrom_size_path_doc,
            mc_contexts_doc=mc_contexts_doc, cov_cutoff_doc=cov_cutoff_doc)
def allc_to_region_count(allc_path: str,
                         output_prefix: str,
                         chrom_size_path: str,
                         mc_contexts: List[str],
                         split_strand: bool = True,
                         region_bed_paths: str = None,
                         region_bed_names: str = None,
                         bin_sizes: List[int] = None,
                         cov_cutoff: int = 9999,
                         save_zero_cov: bool = True,
                         remove_tmp: bool = True):
    """\
    Calculate mC and cov at regional level. Region can be provided in 2 forms:
    1. BED file, provided by region_bed_paths, containing arbitrary regions and use bedtools map to calculate;
    2. Fix-size non-overlap genome bins, provided by bin_sizes.
    Form 2 is much faster to calculate than form 1.
    The output file is in 6-column bed-like format: chrom start end region_uid mc cov

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    output_prefix
        Path prefix of the output region count file.
    chrom_size_path
        {chrom_size_path_doc}
    mc_contexts
        {mc_contexts_doc}
    split_strand
        If true, Watson (+) and Crick (-) strands will be count separately
    region_bed_paths
        Arbitrary genomic regions can be defined in several BED files to count on.
        Space separated paths to each BED files, the fourth column of BED file should be unique id of the region.
    region_bed_names
        Space separated names for each BED file provided in region_bed_paths.
    bin_sizes
        Fix-size genomic bins can be defined by bin_sizes and chrom_size_path.
        Space separated sizes of genome bins, each size will be count separately.
    cov_cutoff
        {cov_cutoff_doc}
    save_zero_cov
        Whether to save the regions that have 0 cov, only apply to region count but not the chromosome count
    remove_tmp
        Whether to remove the temporary BED file

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
                                cov_cutoff=cov_cutoff)
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
