"""
This file is modified from methylpy https://github.com/yupenghe/methylpy

Original author: Yupeng He
"""

import logging
import shlex
import subprocess

from ._open import open_allc
from .utilities import parse_chrom_size

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def _allc_to_bedgraph(allc_path, out_prefix, chrom_size_file,
                      remove_additional_chrom=False, bin_size=50):
    """
    Simply calculate cov and mc_rate for fixed genome bins. No mC context filter.
    """
    chrom_dict = parse_chrom_size(chrom_size_file)
    cur_chrom = 'TOTALLY_NOT_A_CHROM'
    cur_chrom_end = 0
    bin_start = 0
    bin_end = min(cur_chrom_end, bin_size)
    temp_mc, temp_cov = 0, 0

    out_prefix = out_prefix.rstrip('.')
    out_rate = out_prefix + '.rate.bg'
    out_cov = out_prefix + '.cov.bg'

    with open_allc(allc_path) as allc, \
            open(out_rate, 'w') as rate_handle, \
            open(out_cov, 'w') as cov_handle:
        for line in allc:
            chrom, pos, _, _, mc, cov, *_ = line.split("\t")
            pos = int(pos)
            mc = int(mc)
            cov = int(cov)
            if pos >= bin_end or cur_chrom != chrom:
                # write line
                if temp_cov > 0:
                    mc_level = temp_mc / temp_cov
                    rate_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, mc_level])) + "\n")
                    cov_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, temp_cov])) + "\n")

                # reset_chrom
                if cur_chrom != chrom:
                    cur_chrom = chrom
                    try:
                        cur_chrom_end = chrom_dict[chrom]
                    except KeyError as e:
                        # chrom not in chrom size file
                        if remove_additional_chrom:
                            continue
                        else:
                            raise e

                # reset bin
                temp_mc, temp_cov = mc, cov
                bin_start = pos // bin_size * bin_size
                bin_end = min(cur_chrom_end, bin_start + bin_size)
            else:
                temp_mc += mc
                temp_cov += cov
        # write last piece
        if temp_cov > 0:
            mc_level = temp_mc / temp_cov
            rate_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, mc_level])) + "\n")
            cov_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, temp_cov])) + "\n")
    return out_rate, out_cov


def _bedgraph_to_bigwig(input_file, chrom_size_file, path_to_wigtobigwig, remove_bedgraph=True):
    output_file = input_file.rstrip('.bg') + '.bw'
    cmd = f'{path_to_wigtobigwig}wigToBigWig {input_file} {chrom_size_file} {output_file}'

    subprocess.run(shlex.split(cmd), check=True)
    if remove_bedgraph:
        subprocess.run(shlex.split(f'rm -f {input_file}'), check=True)
    return output_file


def allc_to_bigwig(allc_path,
                   output_prefix,
                   chrom_size_path,
                   mc_contexts,
                   strandness='both',
                   bin_size=50,
                   remove_additional_chrom=False,
                   remove_temp_bedgraph=True,
                   path_to_wigtobigwig=""):
    # TODO add context support and strandness support to this function
    # TODO add CLI in __main__
    # TODO write test
    """
    Generate coverage (cov) and ratio (mc/cov) bigwig track files from 1 ALLC file

    Parameters
    ----------
    allc_path
    output_prefix
    chrom_size_path
    bin_size
    remove_additional_chrom
    remove_temp_bedgraph
    path_to_wigtobigwig

    Returns
    -------

    """
    # test wigToBigWig
    p = subprocess.run(f'{path_to_wigtobigwig}wigToBigWig', stderr=subprocess.PIPE, encoding='utf8')
    if p.returncode != 255:
        raise OSError(f'Try {path_to_wigtobigwig}wigToBigWig, got error {p.stderr}')

    # TODO add mc context and strandness split
    # prepare bedgraph
    out_rate, out_cov = _allc_to_bedgraph(allc_path, output_prefix, chrom_size_path,
                                          remove_additional_chrom=remove_additional_chrom,
                                          bin_size=bin_size)

    # generate bigwig file
    _bedgraph_to_bigwig(out_rate, chrom_size_path, path_to_wigtobigwig, remove_temp_bedgraph)
    _bedgraph_to_bigwig(out_cov, chrom_size_path, path_to_wigtobigwig, remove_temp_bedgraph)
    return
