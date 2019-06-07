import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from subprocess import run
from typing import Union, Tuple, Callable, List

import pandas as pd

from ._doc import *
from ._open import open_allc, open_gz
from .utilities import tabix_allc, parse_mc_pattern, parse_chrom_size, genome_region_chunks


def _merge_cg_strand(in_path, out_path):
    """
    Merge strand after extract context step in extract_allc (and only apply on CG), so no need to check context.
    """
    prev_line = None
    cur_chrom = None

    with open_allc(in_path) as allc, \
            open_allc(out_path, 'w') as out_allc:
        for line in allc:
            cur_line = line.strip('\n').split('\t')
            if cur_line[0] != cur_chrom:
                if prev_line is not None:
                    out_allc.write('\t'.join(prev_line) + '\n')
                prev_line = cur_line
                cur_chrom = cur_line[0]
                continue
            if prev_line is None:
                prev_line = cur_line
                continue
            else:
                # pos should be continuous, strand should be reverse
                if int(prev_line[1]) + 1 == int(cur_line[1]) and prev_line[2] != cur_line[2]:
                    new_line = prev_line[:4] + [str(int(prev_line[4]) + int(cur_line[4])),
                                                str(int(prev_line[5]) + int(cur_line[5])), '1']
                    out_allc.write('\t'.join(new_line) + '\n')
                    prev_line = None
                # otherwise, only write and update prev_line
                else:
                    out_allc.write('\t'.join(prev_line) + '\n')
                    prev_line = cur_line
    return None


def _check_strandness_parameter(strandness) -> str:
    strandness = str(strandness).lower()
    if strandness in {'both', 'b'}:
        return 'Both'
    elif strandness in {'merge', 'm'}:
        # first getting both, deal with strand merge later
        return 'MergeTmp'
    elif strandness in {'split', 's'}:
        return 'Split'
    else:
        raise ValueError(f'Unknown value for strandness: {strandness}')


def _check_out_format_parameter(out_format) -> Tuple[str, Callable[[list], str]]:
    def _extract_allc_format(allc_line_list):
        # keep allc format
        return '\t'.join(allc_line_list)

    def _extract_bedgraph_cov_format(allc_line_list):
        # only chrom, pos, pos, cov
        allc_line_list = [allc_line_list[i] for i in [0, 1, 1, 5]]
        return '\t'.join(allc_line_list) + '\n'

    def _extract_bedgraph_rate_format(allc_line_list):
        # only chrom, pos, pos, mc/cov
        allc_line_list = [allc_line_list[i] for i in [0, 1, 1]] + \
                         [f'{(int(allc_line_list[4]) / int(allc_line_list[5])):.3f}']
        return '\t'.join(allc_line_list) + '\n'

    def _extract_bed5_format(allc_line_list):
        # only chrom, pos, pos, mc, cov
        allc_line_list = [allc_line_list[i] for i in [0, 1, 1, 4, 5]]
        return '\t'.join(allc_line_list) + '\n'

    out_format = str(out_format).lower()
    if out_format == 'allc':
        return 'allc.tsv.gz', _extract_allc_format
    elif out_format == 'bed5':
        return 'bed5.bed.gz', _extract_bed5_format
    elif out_format == 'bg-cov':
        return 'cov.bg.gz', _extract_bedgraph_cov_format
    elif out_format == 'bg-rate':
        return 'rate.bg.gz', _extract_bedgraph_rate_format
    else:
        raise ValueError(f'Unknown value for out_format: {out_format}')


def _merge_gz_files(file_list, output_path):
    with open_gz(output_path, 'w') as out_f:
        for file_path in file_list:
            with open_gz(file_path) as f:
                out_f.write(f.read())
            subprocess.run(['rm', '-f', file_path])
    return output_path


def _extract_allc_parallel(allc_path, output_prefix, mc_contexts, strandness, output_format,
                           chrom_size_path, cov_cutoff, cpu, chunk_size=100000000, tabix=True):
    output_prefix = output_prefix.rstrip('.')
    regions = genome_region_chunks(chrom_size_path=chrom_size_path,
                                   bin_length=chunk_size,
                                   combine_small=True)
    future_dict = {}
    with ProcessPoolExecutor(cpu) as executor:
        for chunk_id, region in enumerate(regions):
            future = executor(extract_allc,
                              allc_path=allc_path,
                              output_prefix=output_prefix + f'.{chunk_id}.',
                              mc_contexts=mc_contexts,
                              strandness=strandness,
                              output_format=output_format,
                              chrom_size_path=chrom_size_path,
                              region=region,
                              cov_cutoff=cov_cutoff,
                              cpu=1,
                              tabix=False)
            future_dict[future] = chunk_id

        output_records = []
        for future in as_completed(future_dict):
            output_path_list = future.result()
            output_records += output_path_list

        # agg chunk_output
        records = []
        for path in output_records:
            chunk_id, *full_suffix = path.lstrip(output_prefix).strip('.').split('.')
            full_suffix = '.'.join(full_suffix)
            records.append([path, chunk_id, full_suffix])
        total_output_df = pd.DataFrame(records, columns=['path', 'chunk_id', 'full_suffix'])

        real_out_paths = []
        need_tabix = []
        for suffix, sub_df in total_output_df.groupby('full_suffix'):
            ordered_index = sub_df['chunk_id'].astype(int).sort_values().index
            ordered_file_list = sub_df.loc[ordered_index, 'path'].tolist()
            real_file_path = f'{output_prefix}.{suffix}'
            real_out_paths.append(real_file_path)
            if tabix and 'allc' in suffix:
                need_tabix.append(real_out_paths)
            _merge_gz_files(ordered_file_list, real_file_path)

        if tabix:
            with ProcessPoolExecutor(cpu) as tabix_executor:
                for path in need_tabix:
                    tabix_executor.submit(tabix_allc, path)
    return real_out_paths


@doc_params(allc_path_doc=allc_path_doc,
            mc_contexts_doc=mc_contexts_doc,
            cov_cutoff_doc=cov_cutoff_doc,
            chrom_size_path_doc=chrom_size_path_doc)
def extract_allc(allc_path: str,
                 output_prefix: str,
                 mc_contexts: Union[str, list],
                 chrom_size_path: str,
                 strandness: str = 'both',
                 output_format: str = 'allc',
                 region: str = None,
                 cov_cutoff: int = 9999,
                 tabix: bool = True,
                 cpu=1) -> List[str]:
    """\
    Extract information (strand, context) from 1 ALLC file. Save to several different format.

    Parameters
    ----------
    allc_path
        {allc_path_doc}
    output_prefix
        Path prefix of the output ALLC file.
    mc_contexts
        {mc_contexts_doc}
    strandness
        What to do with strand information, possible values are:
        1. both: save +/- strand together in one file;
        2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-);
        3. merge: This will only merge the count on adjacent CpG in +/- strands, only work for CpG like context.
                  For non-CG context, its the same as both.
    output_format
        Output format of extracted information, possible values are:
        1. allc: keep the allc format
        2. bed5: 5-column bed format, chrom, pos, pos, mc, cov
        3. bg-cov: bedgraph format, chrom, pos, pos, cov
        4. bg-rate: bedgraph format, chrom, pos, pos, mc/cov
    chrom_size_path
        {chrom_size_path_doc}
        If chrom_size_path provided, will use it to extract ALLC with chrom order,
        but if region provided, will ignore this.
    region
        Only extract records from certain genome region(s) via tabix, multiple region can be provided in tabix form.
    cov_cutoff
        {cov_cutoff_doc}
    tabix
        Whether to generate tabix if format is ALLC, only set this to False from _extract_allc_parallel
    cpu
        Number of cores to use

    Returns
    -------
    A list of output file paths, not include index files.
    """
    # TODO write test
    # determine parallel or not
    parallel_chunk_size = 100000000
    if cpu > 1 and region is None:
        return _extract_allc_parallel(allc_path=allc_path,
                                      output_prefix=output_prefix,
                                      mc_contexts=mc_contexts,
                                      strandness=strandness,
                                      output_format=output_format,
                                      chrom_size_path=chrom_size_path,
                                      cov_cutoff=cov_cutoff,
                                      cpu=cpu,
                                      chunk_size=parallel_chunk_size,
                                      tabix=tabix)

    # determine region
    if region is None:
        if chrom_size_path is not None:
            chrom_dict = parse_chrom_size(chrom_size_path)
            region = ' '.join(chrom_dict.keys())

    # prepare params
    output_prefix = output_prefix.rstrip('.')
    if isinstance(mc_contexts, str):
        mc_contexts = mc_contexts.split(' ')
    mc_contexts = list(set(mc_contexts))
    strandness = _check_strandness_parameter(strandness)
    out_suffix, line_func = _check_out_format_parameter(output_format)

    # because mc_contexts can overlap (e.g. CHN, CAN)
    # each context may associate to multiple handle
    context_handle = defaultdict(list)
    handle_collect = []
    output_path_collect = []
    for mc_context in mc_contexts:
        parsed_context_set = parse_mc_pattern(mc_context)
        if strandness == 'Split':
            file_path = output_prefix + f'.{mc_context}-Watson.{out_suffix}'
            output_path_collect.append(file_path)
            w_handle = open_allc(file_path, 'w')
            handle_collect.append(w_handle)

            file_path = output_prefix + f'.{mc_context}-Crick.{out_suffix}'
            output_path_collect.append(file_path)
            c_handle = open_allc(file_path, 'w')
            handle_collect.append(c_handle)
            for mc_pattern in parsed_context_set:
                # handle for Watson/+ strand
                context_handle[(mc_pattern, '+')].append(w_handle)
                # handle for Crick/- strand
                context_handle[(mc_pattern, '-')].append(c_handle)
        else:
            # handle for both strand
            file_path = output_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
            if strandness == 'MergeTmp':
                output_path_collect.append(output_prefix + f'.{mc_context}-Merge.{out_suffix}')
            else:
                output_path_collect.append(file_path)
            _handle = open_allc(file_path, 'w')
            handle_collect.append(_handle)
            for mc_pattern in parsed_context_set:
                context_handle[mc_pattern].append(_handle)

    # split file first
    # strandness function
    with open_allc(allc_path, region=region) as allc:
        if strandness == 'Split':
            for line in allc:
                cur_line = line.split('\t')
                if int(cur_line[5]) > cov_cutoff:
                    continue
                try:
                    # key is (context, strand)
                    [h.write(line_func(cur_line)) for h in context_handle[(cur_line[3], cur_line[2])]]
                except KeyError:
                    continue
        else:
            for line in allc:
                cur_line = line.split('\t')
                if int(cur_line[5]) > cov_cutoff:
                    continue
                try:
                    # key is context
                    [h.write(line_func(cur_line)) for h in context_handle[cur_line[3]]]
                except KeyError:
                    continue
    for handle in handle_collect:
        handle.close()

    if 'allc' in out_suffix:
        for mc_context in mc_contexts:
            # tabix ALLC file
            if strandness == 'Split':
                in_path = output_prefix + f'.{mc_context}-Watson.{out_suffix}'
                if tabix:
                    tabix_allc(in_path)
                in_path = output_prefix + f'.{mc_context}-Crick.{out_suffix}'
                if tabix:
                    tabix_allc(in_path)
            elif strandness == 'MergeTmp':
                in_path = output_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
                if 'CG' in mc_context:
                    out_path = output_prefix + f'.{mc_context}-Merge.{out_suffix}'
                    _merge_cg_strand(in_path, out_path)
                    run(['rm', '-f', in_path], check=True)
                else:
                    # for non-CG, there is no need to merge strand
                    out_path = output_prefix + f'.{mc_context}-Both.{out_suffix}'
                    run(['mv', in_path, out_path], check=True)
                if tabix:
                    tabix_allc(out_path)
            else:
                in_path = output_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
                if tabix:
                    tabix_allc(in_path)
    return output_path_collect
