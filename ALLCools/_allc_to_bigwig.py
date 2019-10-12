"""
This file is modified from methylpy https://github.com/yupenghe/methylpy
"""

import logging
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

from ._doc import *
from ._extract_allc import extract_allc
from ._open import open_allc
from .utilities import parse_chrom_size

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def _allc_to_bedgraph(allc_path, out_prefix, chrom_size_path,
                      remove_additional_chrom=False, bin_size=50):
    """
    Simply calculate cov and mc_rate for fixed genome bins. No mC context filter.
    """

    chrom_dict = parse_chrom_size(chrom_size_path)
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
                # reset_chrom
                if cur_chrom != chrom:
                    try:
                        this_chrom_end = chrom_dict[chrom]
                    except KeyError as e:
                        # chrom not in chrom size file
                        if remove_additional_chrom:
                            continue
                        else:
                            raise e
                    cur_chrom_end = this_chrom_end

                # write line after confirm the chrom
                if temp_cov > 0:
                    mc_level = temp_mc / temp_cov
                    rate_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, mc_level])) + "\n")
                    cov_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, temp_cov])) + "\n")
                cur_chrom = chrom  # only update chrom after wrote the cur line, cause this may belong to last chrom

                # reset bin
                _bin_start = pos // bin_size * bin_size
                if _bin_start == bin_start:
                    # this only happens when pos == cur_chrom_end,
                    # we don't want the last base, just ignore it, usually happens at chrM last base...
                    temp_mc, temp_cov = 0, 0
                else:
                    temp_mc, temp_cov = mc, cov
                    bin_start = _bin_start
                bin_end = min(cur_chrom_end, bin_start + bin_size)
            else:
                temp_mc += mc
                temp_cov += cov

        # write last piece if there is anything
        if temp_cov > 0:
            mc_level = temp_mc / temp_cov
            rate_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, mc_level])) + "\n")
            cov_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end, temp_cov])) + "\n")

    print(f'Finish generate bedgraph for {allc_path}')
    return out_rate, out_cov


def _bedgraph_to_bigwig(input_file, chrom_size_path, path_to_wigtobigwig, remove_bedgraph=True):
    output_file = input_file.rstrip('.bg') + '.bw'
    cmd = f'{path_to_wigtobigwig}wigToBigWig {input_file} {chrom_size_path} {output_file}'

    subprocess.run(shlex.split(cmd), check=True)
    if remove_bedgraph:
        subprocess.run(shlex.split(f'rm -f {input_file}'), check=True)
    print(f'Finish generate bigwig for {input_file}')
    return output_file


@doc_params(allc_path_doc=allc_path_doc,
            chrom_size_path_doc=chrom_size_path_doc,
            mc_contexts_doc=mc_contexts_doc,
            split_strand_doc=split_strand_doc,
            remove_additional_chrom_doc=remove_additional_chrom_doc,
            region_doc=region_doc,
            cov_cutoff_doc=cov_cutoff_doc)
def allc_to_bigwig(allc_path,
                   output_prefix,
                   chrom_size_path,
                   mc_contexts,
                   split_strand=False,
                   bin_size=50,
                   remove_additional_chrom=False,
                   region=None,
                   cov_cutoff=9999,
                   path_to_wigtobigwig="",
                   remove_temp=True,
                   cpu=1):
    """\
    Generate bigwig file(s) from 1 ALLC file.

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    output_prefix
        Output prefix of the bigwig file(s)
    chrom_size_path
        {chrom_size_path_doc}
    mc_contexts
        {mc_contexts_doc}
    split_strand
        {split_strand_doc}
    bin_size
        Minimum bin size of bigwig file
    remove_additional_chrom
        {remove_additional_chrom_doc}
    region
        {region_doc}
    cov_cutoff
        {cov_cutoff_doc}
    path_to_wigtobigwig
        Path to wigtobigwig to allow allcools to find it
    remove_temp
        debug parameter, whether to remove the temp file or not
    cpu
        Number of cores to use

    Returns
    -------

    """
    # TODO write test
    if bin_size is None:
        bin_size = 50

    # test wigToBigWig
    p = subprocess.run(f'{path_to_wigtobigwig}wigToBigWig', stderr=subprocess.PIPE, encoding='utf8')
    if p.returncode != 255:  # somehow, the p.returncode is 255 when set correctly...
        raise OSError(f'Try {path_to_wigtobigwig}wigToBigWig, got error {p.stderr}')

    strandness = 'split' if split_strand else 'both'

    # prepare bedgraph
    extracted_allc_path_dict = extract_allc(allc_path=allc_path,
                                            output_prefix=output_prefix,
                                            mc_contexts=mc_contexts,
                                            chrom_size_path=chrom_size_path,
                                            strandness=strandness,
                                            output_format='allc',
                                            region=region,
                                            cov_cutoff=cov_cutoff,
                                            tabix=False,
                                            cpu=cpu)
    extracted_allc_paths = list(extracted_allc_path_dict.values())

    with ProcessPoolExecutor(cpu) as executor:
        # generate bigwig file
        allc_future_dict = {}
        for path in extracted_allc_paths:
            future = executor.submit(_allc_to_bedgraph,
                                     allc_path=path,
                                     out_prefix=str(path).rstrip('.allc.tsv.gz'),  # use extracted ALLC prefix
                                     chrom_size_path=chrom_size_path,
                                     remove_additional_chrom=remove_additional_chrom,
                                     bin_size=bin_size)
            allc_future_dict[future] = path

        for future in as_completed(allc_future_dict):
            output_paths = future.result()
            for path in output_paths:
                executor.submit(_bedgraph_to_bigwig,
                                input_file=path,
                                chrom_size_path=chrom_size_path,
                                path_to_wigtobigwig=path_to_wigtobigwig,
                                remove_bedgraph=remove_temp)

    if remove_temp:
        subprocess.run(['rm', '-f'] + extracted_allc_paths, check=True)
    return
