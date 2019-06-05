import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import scipy.sparse as ss
from anndata import AnnData

from .utilities import parse_file_paths, \
    parse_chrom_size, \
    chrom_dict_to_id_index, \
    get_bin_id, \
    generate_chrom_bin_bed_dataframe


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
    ss.save_npz(mc_path, mc_matrix, compressed=compression)

    cov_matrix = cov_matrix.tocsr()
    ss.save_npz(cov_path, cov_matrix, compressed=compression)
    return mc_path, cov_path


def _csr_matrix_to_anndata(matrix_paths, obs_names, chrom_size_path, bin_size, output_path,
                           compression=None, compression_opts=None):
    if pathlib.Path(output_path).exists():
        return output_path

    var_df = generate_chrom_bin_bed_dataframe(chrom_size_path, bin_size)
    total_matrix = ss.vstack([ss.load_npz(path) for path in matrix_paths])

    adata = AnnData(X=total_matrix,
                    obs=pd.DataFrame([], index=obs_names),
                    var=var_df[['chrom']],
                    uns=dict(bin_size=bin_size,
                             chrom_size_path=chrom_size_path))
    adata.write(output_path,
                compression=compression,
                compression_opts=compression_opts)
    return output_path


def aggregate_region_count_to_paired_anndata(count_tables, output_prefix, chrom_size_path, bin_size,
                                             compression='gzip', file_uids=None, max_obj=3072, cpu=3):
    # TODO write test
    # this should only deal with a simple case, aggregate 2 sample*feature 2-D matrix, one for mc, one for cov,
    # output to full or sparse format

    mini_batches = min(max_obj, 24)
    output_prefix = output_prefix.rstrip('.')

    # merge bin count table in mini-batches, save each set into npz file
    count_tables = parse_file_paths(count_tables)
    if file_uids is not None:
        file_uids = parse_file_paths(file_uids)
    else:
        file_uids = [pathlib.Path(i).name.rstrip('.sparse.bed.gz') for i in count_tables]
    if len(file_uids) != len(count_tables):
        raise ValueError('Length of file_uids do not match length of count_tables.')

    mc_path_list = []
    cov_path_list = []
    id_chunk_list = []
    future_dict = {}
    worker = cpu // 3
    with ProcessPoolExecutor(worker) as executor:
        for chunk_id, i in enumerate(range(0, len(count_tables), mini_batches)):
            path_chunk = count_tables[i:i + mini_batches]
            cur_id_chunk = file_uids[i:i + mini_batches]
            future = executor.submit(_bin_count_table_to_csr_npz,
                                     bin_count_tables=path_chunk,
                                     bin_size=bin_size,
                                     chrom_size_path=chrom_size_path,
                                     output_prefix=output_prefix + f'.{chunk_id}',
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
    cur_ids = []
    cur_mc_paths = []
    cur_cov_paths = []
    if len(count_tables) > max_obj:
        adata_chunk_id = '.chunk-0'
    else:
        adata_chunk_id = ''
    future_dict = {}
    worker = (cpu - 1) // 6
    with ProcessPoolExecutor(worker) as executor:
        for chunk_id, cur_id_chunk in enumerate(id_chunk_list):
            if len(cur_ids) + len(cur_id_chunk) > max_obj:
                # mc
                executor.submit(_csr_matrix_to_anndata,
                                matrix_paths=cur_mc_paths,
                                obs_names=cur_ids,
                                chrom_size_path=chrom_size_path,
                                bin_size=bin_size,
                                output_path=output_prefix + f'{adata_chunk_id}.mc.h5ad',
                                compression=compression,
                                compression_opts=None)
                # cov
                executor.submit(_csr_matrix_to_anndata,
                                matrix_paths=cur_cov_paths,
                                obs_names=cur_ids,
                                chrom_size_path=chrom_size_path,
                                bin_size=bin_size,
                                output_path=output_prefix + f'{adata_chunk_id}.cov.h5ad',
                                compression=compression,
                                compression_opts=None)
                cur_ids = []
                cur_mc_paths = []
                cur_cov_paths = []
                if adata_chunk_id != '':
                    adata_chunk_id = '.chunk-' + str(int(adata_chunk_id[7:]) + 1)
            else:
                cur_ids += cur_id_chunk
                cur_mc_paths.append(mc_path_list[chunk_id])
                cur_cov_paths.append(cov_path_list[chunk_id])
        if len(cur_ids) != 0:
            # mc
            executor.submit(_csr_matrix_to_anndata,
                            matrix_paths=cur_mc_paths,
                            obs_names=cur_ids,
                            chrom_size_path=chrom_size_path,
                            bin_size=bin_size,
                            output_path=output_prefix + f'{adata_chunk_id}.mc.h5ad',
                            compression=compression,
                            compression_opts=None)
            # cov
            executor.submit(_csr_matrix_to_anndata,
                            matrix_paths=cur_cov_paths,
                            obs_names=cur_ids,
                            chrom_size_path=chrom_size_path,
                            bin_size=bin_size,
                            output_path=output_prefix + f'{adata_chunk_id}.cov.h5ad',
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
