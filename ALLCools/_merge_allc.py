"""
Some of the functions are modified from methylpy https://github.com/yupenghe/methylpy

Original Author: Yupeng He
"""

import gc
import logging
import os
import resource
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import psutil
from pysam import TabixFile
from ._doc import *
from ._open import open_allc
from .utilities import parse_chrom_size, genome_region_chunks, parse_file_paths

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

# get the system soft and hard limit of file handle
SOFT, HARD = resource.getrlimit(resource.RLIMIT_NOFILE)
DEFAULT_MAX_ALLC = 150
PROCESS = psutil.Process(os.getpid())


class _ALLC:
    def __init__(self, path, region):
        self.f = TabixFile(path)
        self.f_region = self.f.fetch(region)

    def readline(self):
        return self.f_region.next()

    def close(self):
        self.f.close()


def _increase_soft_fd_limit():
    """
    Increase soft file descriptor limit to hard limit,
    this is the maximum a process can do
    Use this in merge_allc, because for single cell, a lot of file need to be opened

    Some useful discussion
    https://unix.stackexchange.com/questions/36841/why-is-number-of-open-files-limited-in-linux
    https://docs.python.org/3.6/library/resource.html
    https://stackoverflow.com/questions/6774724/why-python-has-limit-for-count-of-file-handles/6776345
    """
    resource.setrlimit(resource.RLIMIT_NOFILE, (HARD, HARD))


def _batch_merge_allc_files_tabix(allc_files, out_file, chrom_size_file, bin_length, cpu=10, binarize=False, snp=False):
    regions = genome_region_chunks(chrom_size_file, bin_length=bin_length)
    log.info(f'Merge ALLC files with {cpu} processes')
    log.info(f'Split genome into {len(regions)} regions, each is {bin_length}bp')
    log.info(f'{len(allc_files)} to merge, the default ALLC file handel in 1 run is {DEFAULT_MAX_ALLC}')
    log.info(f'Process FH soft limit {SOFT}, hard limit {HARD}')

    _increase_soft_fd_limit()
    if len(allc_files) > DEFAULT_MAX_ALLC:
        # deal with too many allc files
        # merge in batches
        allc_fold = len(allc_files) // DEFAULT_MAX_ALLC + 1
        batch_allc = len(allc_files) // allc_fold + 1

        allc_files_batches = []
        out_paths = []
        for batch_id, i in enumerate(range(0, len(allc_files), batch_allc)):
            allc_files_batches.append(allc_files[i:i + batch_allc])
            out_paths.append(out_file + f'batch_{batch_id}.tmp.tsv.gz')
            if batch_id > 0:
                # merge last batch's merged allc into next batch
                allc_files_batches[batch_id].append(out_paths[batch_id - 1])
        out_paths[-1] = out_file
    else:
        allc_files_batches = [allc_files]
        out_paths = [out_file]

    log.info(f'Run merge in {len(allc_files_batches)} batches')
    log.info(' '.join(out_paths))

    if snp:
        merge_func = _merge_allc_files_tabix_with_snp_info
    else:
        merge_func = _merge_allc_files_tabix

    for batch_num, (allc_files, out_file) in enumerate(zip(allc_files_batches, out_paths)):
        log.info(f'Run batch {batch_num}, '
                 f'{len(allc_files)} allc files, '
                 f'output to {out_file}')
        with open_allc(out_file, 'w', threads=3) as out_handle:
            # as_complete don't release, run total regions in sections to prevent too large memory
            parallel_section = cpu
            for i in range(0, len(regions), parallel_section):
                cur_regions = regions[i:min(i + parallel_section, len(regions))]
                log.info(f'Running region from {cur_regions[0]} to {cur_regions[-1]}')
                with ProcessPoolExecutor(max_workers=cpu) as executor:
                    future_merge_result = {executor.submit(merge_func,
                                                           allc_files=allc_files,
                                                           out_file=None,
                                                           chrom_size_file=chrom_size_file,
                                                           query_region=region,
                                                           buffer_line_number=100000,
                                                           binarize=binarize): region_id
                                           for region_id, region in enumerate(cur_regions)}
                    cur_id = 0
                    temp_dict = {}
                    # future may return in any order
                    # save future.result into temp_dict 1st
                    # write data in order by region_id
                    # so the out file is ordered
                    for future in as_completed(future_merge_result):
                        region_id = future_merge_result[future]
                        try:
                            temp_dict[region_id] = future.result()
                        except Exception as exc:
                            log.info('%r generated an exception: %s' % (region_id, exc))
                        else:
                            try:
                                while len(temp_dict) > 0:
                                    data = temp_dict.pop(cur_id)
                                    out_handle.write(data)
                                    log.info(f'write {cur_regions[cur_id]} Cached: {len(temp_dict)}, '
                                             f'Current memory size: {PROCESS.memory_info().rss / (1024 ** 3):.2f}')
                                    cur_id += 1
                            except KeyError:
                                continue
                    # write last pieces of data
                    while len(temp_dict) > 0:
                        if cur_id in temp_dict:
                            data = temp_dict.pop(cur_id)
                            out_handle.write(data)
                            log.info(f'write {regions[cur_id]} Cached: {len(temp_dict)}, '
                                     f'Current memory size: {PROCESS.memory_info().rss / (1024 ** 3):.2f}')
                        cur_id += 1
                gc.collect()
        # after merge, tabix output
        log.info('Tabix output ALLC file')
        subprocess.run(['tabix', '-b', '2', '-e', '2', '-s', '1', out_file],
                       check=True)
        log.info(f'Current memory size: {PROCESS.memory_info().rss / (1024 ** 3):.2f}')
    log.info('Merge finished.')

    # remove all batch allc but the last (final allc)
    for out_file in out_paths[:-1]:  # last file is the final merged allc
        subprocess.run(shlex.split(f'rm -f {out_file} {out_file}.tbi'))
    return


def _merge_allc_files_tabix(allc_files,
                            out_file,
                            chrom_size_file,
                            query_region=None,
                            buffer_line_number=10000,
                            binarize=False):
    # only use bgzip and tabix
    # automatically take care the too many file open issue
    # do merge iteratively if file number exceed limit
    # parallel in chrom_bin level, not chrom level

    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")
    chrom_size_dict = parse_chrom_size(chrom_size_file)
    all_chroms = list(chrom_size_dict.keys())
    if query_region is not None:
        if not isinstance(query_region, str):
            exit("query_region must be str or None")
        region_chroms = set([region.split(':')[0] for region in query_region.split(' ')])
        all_chroms = [chrom for chrom in all_chroms if chrom in region_chroms]
    processing_chrom = all_chroms.pop(0)

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    file_handles = [_ALLC(allc_file,
                          region=query_region)
                    for allc_file in allc_files]
    if out_file is not None:
        out_handle = open_allc(out_file, 'w', threads=3)
    else:
        out_handle = ''

    # merge allc files
    out = ""
    cur_chrom = ['NOT_A_CHROM' for _ in range(len(allc_files))]
    cur_pos = np.array([np.nan for _ in range(len(allc_files))])
    cur_fields = [None for _ in range(len(allc_files))]
    file_reading = np.array([True for _ in range(len(allc_files))])

    # init
    for index, allc_file in enumerate(allc_files):
        try:
            line = file_handles[index].readline()
            fields = line.split("\t")
            cur_chrom[index] = fields[0]
            cur_pos[index] = int(fields[1])
            cur_fields[index] = fields
        except StopIteration:
            # file handle read nothing, the file is empty
            file_reading[index] = False

    active_handle = np.array([True if chrom == processing_chrom else False
                              for chrom in cur_chrom])

    # merge
    line_count = 0
    while file_reading.sum() > 0:
        mc, cov = 0, 0
        genome_info = None
        # select index whose cur_pos is smallest among all active handle
        for index in np.where((cur_pos == np.nanmin(cur_pos[active_handle]))
                              & active_handle)[0]:
            mc += int(cur_fields[index][4])
            cov += int(cur_fields[index][5])
            # TODO fix binarize ALLC problem
            # if binarize:
            #     mc, cov = binary_count(int(mc), int(cov))
            #     mc += mc
            #     cov += cov
            # else:
            #     mc += int(cur_fields[index][4])
            #     cov += int(cur_fields[index][5])
            if genome_info is None:
                genome_info = cur_fields[index][:4]

            # update
            try:
                line = file_handles[index].readline()
                fields = line.split("\t")
            except StopIteration:
                # read to the end of a file
                fields = ['NOT_A_CHROM', 9999999999]
                file_reading[index] = False

            # judge if chrom changed between two lines
            this_chrom = cur_fields[index][0]
            next_chrom = fields[0]
            if next_chrom == this_chrom:
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields
            else:
                # read to next chrom
                # handle became inactive
                active_handle[index] = False
                cur_chrom[index] = next_chrom
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields

                # if all handle became inactive, move processing_chrom to next
                if sum(active_handle) == 0:
                    if len(all_chroms) == 0:
                        break
                    processing_chrom = all_chroms.pop(0)
                    # and re-judge active handle
                    active_handle = np.array([True if chrom == processing_chrom else False
                                              for chrom in cur_chrom])

        # output
        out += '\t'.join(genome_info) + f'\t{mc}\t{cov}\t1\n'
        line_count += 1
        if line_count > buffer_line_number:
            if isinstance(out_handle, str):
                out_handle += out
            else:
                out_handle.write(out)
            line_count = 0
            out = ""
    # the last out
    for file_handle in file_handles:
        file_handle.close()
    if isinstance(out_handle, str):
        out_handle += out
        return out_handle
    else:
        out_handle.write(out)
        out_handle.close()
        return


def _merge_allc_files_tabix_with_snp_info(allc_files,
                                          out_file,
                                          chrom_size_file,
                                          query_region=None,
                                          buffer_line_number=10000,
                                          binarize=False):
    # only use bgzip and tabix
    # automatically take care the too many file open issue
    # do merge iteratively if file number exceed limit
    # parallel in chrom_bin level, not chrom level

    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")
    chrom_size_dict = parse_chrom_size(chrom_size_file)
    all_chroms = list(chrom_size_dict.keys())
    if query_region is not None:
        if not isinstance(query_region, str):
            exit("query_region must be str or None")
        region_chroms = set([region.split(':')[0] for region in query_region.split(' ')])
        all_chroms = [chrom for chrom in all_chroms if chrom in region_chroms]
    processing_chrom = all_chroms.pop(0)

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    file_handles = [_ALLC(allc_file,
                          region=query_region)
                    for allc_file in allc_files]
    if out_file is not None:
        out_handle = open_allc(out_file, 'w')
    else:
        out_handle = ''

    # merge allc files
    out = ""
    cur_chrom = ['NOT_A_CHROM' for _ in range(len(allc_files))]
    cur_pos = np.array([np.nan for _ in range(len(allc_files))])
    cur_fields = [None for _ in range(len(allc_files))]
    file_reading = np.array([True for _ in range(len(allc_files))])

    # init
    for index, allc_file in enumerate(allc_files):
        line = file_handles[index].readline()
        if line:
            fields = line.split("\t")
            cur_chrom[index] = fields[0]
            cur_pos[index] = int(fields[1])
            cur_fields[index] = fields
        else:
            # file handle read nothing, the file is empty
            file_reading[index] = False

    active_handle = np.array([True if chrom == processing_chrom else False
                              for chrom in cur_chrom])

    # merge
    line_count = 0
    while file_reading.sum() > 0:
        mc, cov = 0, 0
        # snp specific
        snp_match, snp_mismatch = [0, 0, 0], [0, 0, 0]
        genome_info = None
        # select index whose cur_pos is smallest among all active handle
        for index in np.where((cur_pos == np.nanmin(cur_pos[active_handle]))
                              & active_handle)[0]:
            mc += int(cur_fields[index][4])
            cov += int(cur_fields[index][5])

            # snp specific
            # col 7 is match, col 8 is mismatch
            snp_match = list(map(sum, zip(snp_match, map(int, cur_fields[index][7].split(',')))))
            snp_mismatch = list(map(sum, zip(snp_mismatch, map(int, cur_fields[index][8].split(',')))))

            # TODO fix binarize ALLC problem
            # if binarize:
            #     mc, cov = binary_count(int(mc), int(cov))
            #     mc += mc
            #     cov += cov
            # else:
            #     mc += int(cur_fields[index][4])
            #     cov += int(cur_fields[index][5])
            if genome_info is None:
                genome_info = cur_fields[index][:4]

            # update
            try:
                line = file_handles[index].readline()
                fields = line.split("\t")
            except StopIteration:
                # read to the end of a file
                fields = ['NOT_A_CHROM', 9999999999]
                file_reading[index] = False

            # judge if chrom changed between two lines
            this_chrom = cur_fields[index][0]
            next_chrom = fields[0]
            if next_chrom == this_chrom:
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields
            else:
                # read to next chrom
                # handle became inactive
                active_handle[index] = False
                cur_chrom[index] = next_chrom
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields

                # if all handle became inactive, move processing_chrom to next
                if sum(active_handle) == 0:
                    if len(all_chroms) == 0:
                        break
                    processing_chrom = all_chroms.pop(0)
                    # and re-judge active handle
                    active_handle = np.array([True if chrom == processing_chrom else False
                                              for chrom in cur_chrom])

        # output
        out += '\t'.join(
            genome_info) + f'\t{mc}\t{cov}\t1\t{",".join(map(str, snp_match))}\t{",".join(map(str, snp_mismatch))}\n'
        line_count += 1
        if line_count > buffer_line_number:
            if isinstance(out_handle, str):
                out_handle += out
            else:
                out_handle.write(out)
            line_count = 0
            out = ""
    # the last out
    for file_handle in file_handles:
        file_handle.close()
    if isinstance(out_handle, str):
        out_handle += out
        return out_handle
    else:
        out_handle.write(out)
        out_handle.close()
        return


@doc_params(allc_paths_doc=allc_paths_doc,
            chrom_size_path_doc=chrom_size_path_doc,
            cpu_basic_doc=cpu_basic_doc,
            binarize_doc=binarize_doc,
            snp_doc=snp_doc)
def merge_allc_files(allc_paths, output_path, chrom_size_path, bin_length=10000000, cpu=10, binarize=False, snp=False):
    """
    Merge N ALLC files into 1 ALLC file.

    Parameters
    ----------
    allc_paths
        {allc_paths_doc}
    output_path
        Path to the output merged ALLC file.
    chrom_size_path
        {chrom_size_path_doc}
    bin_length
        Length of the genome bin in each parallel job, large number means more memory usage.
    cpu
        {cpu_basic_doc}
        The real CPU usage is ~1.5 times than this number,
        due to the sub processes of handling ALLC files using tabix/bgzip.
        Monitor the CPU and Memory usage when running this function.
    binarize
        {binarize_doc}
    snp
        {snp_doc}
    Returns
    -------

    """
    # TODO binarize do not work because when merge batch allc, the previous batch merged allc is not single cell
    # can not treat that as binarize, need to reimplement merge ALLC to support binarize in merge allc
    log.info('Right now binarize is not used, need fix this in merge ALLC fiction, set binarize=False')
    binarize = False

    # TODO write test
    # a list of str, contain all absolute non-soft-link paths
    allc_files: list = parse_file_paths(allc_paths)
    if len(allc_files) < 2:
        raise ValueError('There is less than 2 files after parsing the provided allc_paths.')

    try:
        with open(output_path, 'w'):
            pass
    except IOError:
        log.info("Can't create output_path")

    try:
        for allc_path in allc_files:
            if not os.path.exists(allc_path + '.tbi'):
                raise FileNotFoundError(f'Tabix for {allc_path} not found')
    except FileNotFoundError:
        log.info('Some ALLC file do not have tabix, in order to use this function, '
                 'you need to bgzip compress the ALLC file and use tabix to generate .tbi index. (ALLCools defaults)')

    _batch_merge_allc_files_tabix(allc_files=allc_files,
                                  out_file=output_path,
                                  chrom_size_file=chrom_size_path,
                                  bin_length=bin_length,
                                  cpu=cpu,
                                  binarize=binarize,
                                  snp=snp)
    return
