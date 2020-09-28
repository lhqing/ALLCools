import pathlib
import shlex
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import List

import pandas as pd

from ._doc import *
from ._extract_allc import extract_allc
from ._open import open_gz
from .utilities import \
    check_tbi_chroms, \
    parse_chrom_size, \
    chrom_dict_to_id_index, \
    get_bin_id, \
    _transfer_bin_size, \
    generate_chrom_bin_bed_dataframe


def _bedtools_map(region_bed, site_bed, out_bed, chrom_size_path, save_zero_cov=True):
    """
    Use bedtools map to map site_bed format into any bed file provided.
    """
    cmd = f'bedtools map -a {region_bed} -b {site_bed} -c 4,5 -o sum,sum -null 0 -g {chrom_size_path}'
    print(cmd)
    bed_out = subprocess.Popen(shlex.split(cmd),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               encoding='utf8')
    if out_bed.endswith('gz'):
        opener = partial(open_gz, mode='wt')
    else:
        opener = partial(open, mode='w')

    with opener(out_bed) as out_handle:
        while True:
            line = bed_out.stdout.readline()
            if line == '':
                break
            if save_zero_cov or (not line.endswith('\t0\n')):
                out_handle.write(line)

    return_code = bed_out.poll()
    if return_code is None:
        # still running
        bed_out.terminate()
        time.sleep(0.01)
    return_code = bed_out.poll()
    if return_code is not None and return_code != 0:
        message = bed_out.stderr.read().strip()
        raise OSError(message)

    return out_bed


def _map_to_sparse_chrom_bin(site_bed, out_bed, chrom_size_path,
                             bin_size=500):
    """
    Calculate chromosome bins regional count, output is SPARSE,
    bin_id constructed from chrom_size_path and can be reproduce.
    """
    chrom_dict = parse_chrom_size(chrom_size_path)
    chrom_index_dict = chrom_dict_to_id_index(chrom_dict, bin_size)
    cur_chrom = 'TOTALLY_NOT_A_CHROM'
    cur_chrom_end = 0
    bin_end = min(cur_chrom_end, bin_size)
    bin_start = 0
    bin_id = -1
    temp_mc, temp_cov = -1, -1

    with open_gz(site_bed) as in_handle, \
            open_gz(out_bed, 'w') as out_handle:
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
                bin_id = get_bin_id(cur_chrom, chrom_index_dict, bin_start, bin_size)
            else:
                temp_mc += mc
                temp_cov += cov
        # write last piece
        if temp_cov > 0:
            out_handle.write("\t".join(map(str, [cur_chrom, bin_start, bin_end,
                                                 bin_id, temp_mc, temp_cov])) + "\n")
    return out_bed


@doc_params(allc_path_doc=allc_path_doc,
            chrom_size_path_doc=chrom_size_path_doc,
            mc_contexts_doc=mc_contexts_doc,
            cov_cutoff_doc=cov_cutoff_doc,
            split_strand_doc=split_strand_doc,
            cpu_basic_doc=cpu_basic_doc,
            bin_sizes_doc=bin_sizes_doc,
            region_bed_paths_doc=region_bed_paths_doc,
            region_bed_names_doc=region_bed_names_doc,
            binarize_doc=binarize_doc)
def allc_to_region_count(allc_path: str,
                         output_prefix: str,
                         chrom_size_path: str,
                         mc_contexts: List[str],
                         split_strand: bool = False,
                         region_bed_paths: List[str] = None,
                         region_bed_names: List[str] = None,
                         bin_sizes: List[int] = None,
                         cov_cutoff: int = 9999,
                         save_zero_cov: bool = False,
                         remove_tmp: bool = True,
                         cpu: int = 1,
                         binarize: bool = False):
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
        {split_strand_doc}
    region_bed_paths
        {region_bed_paths_doc}
    region_bed_names
        {region_bed_names_doc}
    bin_sizes
        {bin_sizes_doc}
    cov_cutoff
        {cov_cutoff_doc}
    save_zero_cov
        Whether to save the regions that have 0 cov, only apply to region count but not the chromosome count
    remove_tmp
        Whether to remove the temporary BED file
    cpu
        {cpu_basic_doc}
        This function parallel on region level at the extraction step
        and will generate a bunch of small files if cpu > 1.
        Do not use cpu > 1 for single cell region count.
        For single cell data, parallel on cell level is better.
    binarize
        {binarize_doc}
    Returns
    -------
    """
    # TODO write test
    cpu = int(cpu)

    genome_dict = parse_chrom_size(chrom_size_path)
    if bin_sizes is None and region_bed_paths is None:
        raise ValueError('Either bin_sizes or region_bed_paths should be provided.')

    # check bed file
    # 1. bgzip and tabix
    # 2. order of chrom should be the same as order of chrom_size_path
    if region_bed_paths is not None:
        if (region_bed_names is None) or (len(region_bed_names) != len(region_bed_paths)):
            raise ValueError('Different number of bed names and paths provided')
        for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
            if not check_tbi_chroms(region_bed_path, genome_dict):
                raise ValueError(f'Make sure the bed file {region_bed_path} is:'
                                 f'1. bgzip and tabixed;'
                                 f'2. chromosome order is same as the {chrom_size_path}.')

    # print('Extract ALLC context')
    output_prefix = output_prefix.rstrip('.')
    strandness = 'split' if split_strand else 'both'
    output_paths_dict = extract_allc(allc_path=allc_path,
                                     output_prefix=output_prefix,
                                     mc_contexts=mc_contexts,
                                     strandness=strandness,
                                     output_format='bed5',
                                     chrom_size_path=chrom_size_path,
                                     region=None,
                                     cov_cutoff=cov_cutoff,
                                     cpu=cpu,
                                     binarize=binarize)

    path_dict = {}
    for (mc_context, strandness, _), path in output_paths_dict.items():
        # this is according to extract_allc return format
        info_type = f'{mc_context}-{strandness}'
        path_dict[info_type] = path

    with ProcessPoolExecutor(cpu) as executor:
        futures = []
        if region_bed_paths is not None:
            # print('Map to regions.')
            save_flag = 'full' if save_zero_cov else 'sparse'
            for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
                for info_type, site_bed_path in path_dict.items():
                    future = executor.submit(_bedtools_map,
                                             region_bed=region_bed_path,
                                             site_bed=site_bed_path,
                                             out_bed=output_prefix + f'.{region_name}'
                                             f'_{info_type}.{save_flag}.bed.gz',
                                             chrom_size_path=chrom_size_path,
                                             save_zero_cov=save_zero_cov)
                    futures.append(future)

        if bin_sizes is not None:
            # print('Map to chromosome bins.')
            for bin_size in bin_sizes:
                for info_type, site_bed_path in path_dict.items():
                    future = executor.submit(_map_to_sparse_chrom_bin,
                                             site_bed=site_bed_path,
                                             out_bed=output_prefix + f'.chrom{_transfer_bin_size(bin_size)}'
                                             f'_{info_type}.sparse.bed.gz',
                                             chrom_size_path=chrom_size_path,
                                             bin_size=bin_size)
                    futures.append(future)
        output_collection = []
        for future in as_completed(futures):
            # call all future.result to make sure finished successfully
            output_collection.append(future.result())

    if remove_tmp:
        # print('Remove temporal files.')
        for site_bed_path in path_dict.values():
            subprocess.run(['rm', '-f', site_bed_path])
            subprocess.run(['rm', '-f', site_bed_path + '.tbi'])

    # TODO collect all output path, return a informative dict
    return output_collection


def batch_allc_to_region_count(allc_series,
                               output_dir,
                               chrom_size_path,
                               mc_contexts,
                               split_strand,
                               bin_sizes=None,
                               region_bed_paths=None,
                               region_bed_names=None,
                               cov_cutoff=9999,
                               cpu=5,
                               binarize=False):
    output_dir = pathlib.Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    # dump region_bed_path to output_dir for future records
    if region_bed_paths is not None:
        for region_bed_name, region_bed_path in zip(region_bed_names, region_bed_paths):
            bed_df = pd.read_csv(region_bed_path, header=None, index_col=3, sep='\t')
            bed_df.columns = ['chrom', 'start', 'end']
            bed_df.index.name = region_bed_name
            bed_df['int_id'] = list(range(0, bed_df.shape[0]))
            bed_df['int_id'].to_hdf(output_dir / f'REGION_ID_{region_bed_name}.hdf', key='data')
            bed_df.iloc[:, :3].to_hdf(output_dir / f'REGION_BED_{region_bed_name}.hdf', key='data')

    if bin_sizes is not None:
        for bin_size in bin_sizes:
            bin_size_chr = _transfer_bin_size(bin_size)
            region_name = f'chrom{bin_size_chr}'
            bed_df = generate_chrom_bin_bed_dataframe(chrom_size_path=chrom_size_path,
                                                      window_size=bin_size,
                                                      step_size=bin_size)
            bed_df['int_id'] = list(range(0, bed_df.shape[0]))
            bed_df['int_id'].to_hdf(output_dir / f'REGION_ID_{region_name}.hdf', key='data')
            bed_df.iloc[:, :3].to_hdf(output_dir / f'REGION_BED_chrom{bin_size_chr}.hdf', key='data')

    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for cell_id, allc_path in allc_series.iteritems():
            future = executor.submit(allc_to_region_count,
                                     allc_path=allc_path,
                                     output_prefix=str(output_dir / cell_id),
                                     chrom_size_path=chrom_size_path,
                                     mc_contexts=mc_contexts,
                                     split_strand=split_strand,
                                     region_bed_paths=region_bed_paths,
                                     region_bed_names=region_bed_names,
                                     bin_sizes=bin_sizes,
                                     cov_cutoff=cov_cutoff,
                                     save_zero_cov=False,
                                     remove_tmp=True,
                                     cpu=1,
                                     binarize=binarize)
            futures[future] = cell_id

        records = {}
        finish_count = 0
        for future in as_completed(futures):
            cell_id = futures[future]
            try:
                output_path_collect = future.result()
                records[cell_id] = output_path_collect
            except Exception as e:
                print(f'{cell_id} raised an error!')
                raise e
            finish_count += 1
            if finish_count % 100 == 0:
                print(f'{finish_count} allc-to-region-count finished.')

    path_records = []
    for file_id, out_paths in records.items():
        for path in out_paths:
            file_name = pathlib.Path(path).name
            *_, region_mc_type_strand, _, _, _ = file_name.split('.')
            region_name, mc_type, strandness = region_mc_type_strand.replace('_', '-').split('-')
            path_dict = {
                'region_name': region_name,
                'mc_type': mc_type,
                'strandness': strandness,
                'file_id': file_id,
                'file_path': path
            }
            path_records.append(path_dict)
    path_df = pd.DataFrame(path_records)
    path_df.to_hdf(output_dir / 'REGION_COUNT_SUMMARY.hdf', key='data')
    return
