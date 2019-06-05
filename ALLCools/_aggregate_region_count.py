import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import scipy.sparse as ss

from anndata import h5py
from ._mcad_io import csr_matrix_to_mcad
from .utilities import parse_file_paths, \
    parse_chrom_size, \
    chrom_dict_to_id_index, \
    get_bin_id


def _bin_count_table_to_csr_npz(bin_count_tables,
                                bin_size,
                                chrom_size_path,
                                output_prefix,
                                compression=True):
    output_prefix = output_prefix.rstrip('.')
    mc_path = output_prefix + '.mc.npz'
    cov_path = output_prefix + '.cov.npz'
    if pathlib.Path(mc_path).exists() and pathlib.Path(cov_path).exists():
        return mc_path, cov_path

    count_table_path_list = parse_file_paths(bin_count_tables)
    chrom_size_dict = parse_chrom_size(chrom_size_path)
    chrom_bin_index_dict = chrom_dict_to_id_index(chrom_size_dict, bin_size)

    # determine matrix shape
    n_obj = len(count_table_path_list)
    last_bin = get_bin_id(chrom=list(chrom_size_dict.keys())[-1],
                          chrom_index_dict=chrom_bin_index_dict,
                          bin_start=list(chrom_size_dict.values())[-1],
                          bin_size=bin_size)
    n_feature = last_bin + 1

    # init matrix
    # hint by scipy, lil_matrix is the only suitable format in construction step
    mc_matrix = ss.lil_matrix((n_obj, n_feature), dtype=np.uint32)
    cov_matrix = ss.lil_matrix((n_obj, n_feature), dtype=np.uint32)

    # read each count table and add to matrix
    for obj_idx, file_path in enumerate(count_table_path_list):
        data = pd.read_csv(file_path,
                           usecols=[3, 4, 5],
                           sep='\t',
                           header=None,
                           names=['bin_id', 'mc', 'cov'],
                           dtype={'bin_id': np.uint64,
                                  'mc': np.uint32,
                                  'cov': np.uint32},
                           index_col=0)
        mc_matrix[obj_idx, data.index.values] = data['mc'].values
        cov_matrix[obj_idx, data.index.values] = data['cov'].values
    mc_matrix = mc_matrix.tocsr()
    ss.save_npz(mc_path + '.temp.npz', mc_matrix, compressed=compression)
    subprocess.run(['mv', mc_path + '.temp.npz', mc_path], check=True)  # remove temp to make sure file is intact
    cov_matrix = cov_matrix.tocsr()
    ss.save_npz(cov_path + '.temp.npz', cov_matrix, compressed=compression)
    subprocess.run(['mv', cov_path + '.temp.npz', cov_path], check=True)
    return mc_path, cov_path


def aggregate_region_count_to_mcad(count_tables,
                                   output_path,
                                   chrom_size_path,
                                   bin_size,
                                   mc_type: str,
                                   strandness: str,
                                   compression='gzip',
                                   file_uids: list = None,
                                   max_obj=3072,
                                   mini_batches=48,
                                   cpu=5):
    # TODO write CLI, doc and test
    # mcad structure
    # /{mc_type}-{strandness}/  # information type
    #   0/  # cell chunk id
    #     mc/  # count type
    #       X/  # anndata-like structure
    #       obs/
    #       var/
    #       ...
    #     cov/
    #       X/
    #       obx/
    #       var/
    #       ...
    #   1/  # next chunk (optional)

    # deal with user input
    strandness = strandness.capitalize()
    if not strandness in ['Both', 'Split', 'Watson', 'Crick']:
        raise ValueError(f'Unknown strandness value: {strandness}')
    information_type = f'{mc_type.upper()}-{strandness}'
    mini_batches = min(max_obj, mini_batches)
    if pathlib.Path(output_path).exists():
        raise FileExistsError(f'{output_path} already exists')
    if not output_path.endswith('.mcad'):
        raise ValueError('output_path should end with .mcad')

    # merge bin count table in mini-batches, save each set into 2 npz file, one for mc, one for cov
    count_tables = parse_file_paths(count_tables)
    if file_uids is None:
        file_uids = [pathlib.Path(i).name.rstrip('.sparse.bed.gz') for i in count_tables]
    if len(set(file_uids)) != len(file_uids):
        raise ValueError('file_uids is not unique.')
    if len(file_uids) != len(count_tables):
        raise ValueError('Length of file_uids do not match length of count_tables.')

    mc_path_list = []
    cov_path_list = []
    id_chunk_list = []
    worker = (cpu - 1) // 2  # because numpy and pandas use multi-threads already
    future_dict = {}
    with ProcessPoolExecutor(worker) as executor:
        for chunk_id, i in enumerate(range(0, len(count_tables), mini_batches)):
            path_chunk = count_tables[i:i + mini_batches]
            cur_id_chunk = file_uids[i:i + mini_batches]
            future = executor.submit(_bin_count_table_to_csr_npz,
                                     bin_count_tables=path_chunk,
                                     bin_size=bin_size,
                                     chrom_size_path=chrom_size_path,
                                     output_prefix=output_path + f'.{chunk_id}',
                                     compression=True)
            future_dict[future] = cur_id_chunk

        for future in as_completed(future_dict):
            cur_id_chunk = future_dict[future]
            try:
                mc_path, cov_path = future.result()
                mc_path_list.append(mc_path)
                cov_path_list.append(cov_path)
                id_chunk_list.append(cur_id_chunk)
            except Exception as e:
                print(f'Got error when running chunk id list: {cur_id_chunk}')
                raise e

    # merge all npz file into AnnData
    # init
    cur_ids = []
    cur_mc_paths = []
    cur_cov_paths = []
    adata_chunk_id = 0
    # writing to hdf5 is not able to parallel
    with h5py.File(output_path, 'a') as f:
        for chunk_id, cur_id_chunk in enumerate(id_chunk_list):
            if len(cur_ids) + len(cur_id_chunk) > max_obj:
                # mc
                csr_matrix_to_mcad(key_base=f'{information_type}/{adata_chunk_id}/mc',
                                   matrix_paths=cur_mc_paths,
                                   obs_names=cur_ids,
                                   chrom_size_path=chrom_size_path,
                                   bin_size=bin_size,
                                   f=f,
                                   compression=compression,
                                   compression_opts=None)
                # cov
                csr_matrix_to_mcad(key_base=f'{information_type}/{adata_chunk_id}/cov',
                                   matrix_paths=cur_cov_paths,
                                   obs_names=cur_ids,
                                   chrom_size_path=chrom_size_path,
                                   bin_size=bin_size,
                                   f=f,
                                   compression=compression,
                                   compression_opts=None)
                cur_ids = []
                cur_mc_paths = []
                cur_cov_paths = []
                adata_chunk_id += 1
            else:
                cur_ids += cur_id_chunk
                cur_mc_paths.append(mc_path_list[chunk_id])
                cur_cov_paths.append(cov_path_list[chunk_id])
        if len(cur_ids) != 0:
            # mc
            csr_matrix_to_mcad(key_base=f'{information_type}/{adata_chunk_id}/mc',
                               matrix_paths=cur_mc_paths,
                               obs_names=cur_ids,
                               chrom_size_path=chrom_size_path,
                               bin_size=bin_size,
                               f=f,
                               compression=compression,
                               compression_opts=None)
            # cov
            csr_matrix_to_mcad(key_base=f'{information_type}/{adata_chunk_id}/cov',
                               matrix_paths=cur_cov_paths,
                               obs_names=cur_ids,
                               chrom_size_path=chrom_size_path,
                               bin_size=bin_size,
                               f=f,
                               compression=compression,
                               compression_opts=None)
        # remove temp file until everything finished
    subprocess.run(['rm', '-f'] + mc_path_list + cov_path_list)
    return


def aggregate_region_count_to_mcds(count_table_dataframe):
    # TODO write prepare mcds
    # TODO write test
    # MCDS dimension: cell, feature, count_type, mc_type, strand
    # count_table_dataframe: at least cell, feature, file_path, optional mc_type, strand
    return
