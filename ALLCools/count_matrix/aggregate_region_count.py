import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import scipy.sparse as ss
import xarray as xr
from anndata import AnnData

from .._allc_to_region_count import allc_to_region_count
from ..utilities import parse_file_paths, \
    parse_chrom_size, \
    chrom_dict_to_id_index, \
    get_bin_id, \
    generate_chrom_bin_bed_dataframe, \
    _transfer_bin_size


def _batch_allc_to_region_count(allc_table,
                                output_dir,
                                chrom_size_path,
                                mc_contexts,
                                split_strand,
                                bin_sizes,
                                region_bed_paths=None,
                                region_bed_names=None,
                                cpu=5):
    if isinstance(allc_table, str):
        allc_series = pd.read_csv(allc_table, header=None, index_col=0, squeeze=True, sep='\t')
        if not isinstance(allc_series, pd.Series):
            raise ValueError('allc_table malformed, should only have 2 columns, 1. file_uid, 2. file_path')
    else:
        allc_series = allc_table
    if allc_series.index.duplicated().sum() != 0:
        raise ValueError('allc_table file uid have duplicates (1st column)')

    output_dir = pathlib.Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    # dump region_bed_path to output_dir for future records
    if region_bed_paths is not None:
        for region_bed_name, region_bed_path in zip(region_bed_names, region_bed_paths):
            bed_df = pd.read_csv(region_bed_path, header=None, index_col=3, sep='\t')
            bed_df['int_id'] = list(range(0, bed_df.shape[0]))
            bed_df['int_id'].to_msgpack(output_dir / f'REGION_ID_{region_bed_name}.msg')

    if bin_sizes is not None:
        for bin_size in bin_sizes:
            bin_size_chr = _transfer_bin_size(bin_size)
            region_name = f'chrom{bin_size_chr}'
            bed_df = generate_chrom_bin_bed_dataframe(chrom_size_path=chrom_size_path,
                                                      window_size=bin_size,
                                                      step_size=bin_size)
            bed_df['int_id'] = list(range(0, bed_df.shape[0]))
            bed_df['int_id'].to_msgpack(output_dir / f'REGION_ID_{region_name}.msg')

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
                                     cov_cutoff=2,
                                     save_zero_cov=False,
                                     remove_tmp=True,
                                     cpu=1)
            futures[future] = cell_id

        records = {}
        for future in as_completed(futures):
            cell_id = futures[future]
            try:
                output_path_collect = future.result()
                records[cell_id] = output_path_collect
            except Exception as e:
                print(f'{cell_id} raised an error!')
                raise e

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
    path_df.to_msgpack(output_dir / 'REGION_COUNT_SUMMARY.msg')
    return


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


def _csr_matrix_to_anndata(matrix_paths, output_path, obs_names, chrom_size_path,
                           bin_size, mc_type, count_type, step_size, strandness,
                           compression=None, compression_opts=None):
    if pathlib.Path(output_path).exists():
        return output_path

    var_df = generate_chrom_bin_bed_dataframe(chrom_size_path, bin_size)
    total_matrix = ss.vstack([ss.load_npz(path) for path in matrix_paths])

    adata = AnnData(X=total_matrix,
                    obs=pd.DataFrame([], index=obs_names),
                    var=var_df[['chrom']],
                    uns=dict(bin_size=bin_size,
                             chrom_size_path=chrom_size_path,
                             mc_type=mc_type,
                             count_type=count_type,
                             step_size=step_size,
                             strandness=strandness))
    adata.write(output_path,
                compression=compression,
                compression_opts=compression_opts)
    return output_path


def aggregate_region_count_to_paired_anndata(count_tables, output_prefix, chrom_size_path,
                                             bin_size, mc_type, count_type, strandness,
                                             compression='gzip', file_uids=None, max_obj=3072, cpu=3):
    # TODO write test
    # this should only deal with a simple case, aggregate 2 sample*feature 2-D matrix, one for mc, one for cov,
    # output to full or sparse format
    step_size = bin_size  # which means the region is non-overlap

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
                                mc_type=mc_type,
                                count_type=count_type,
                                step_size=step_size,
                                strandness=strandness,
                                output_path=output_prefix + f'{adata_chunk_id}.mc.h5ad',
                                compression=compression,
                                compression_opts=None)
                # cov
                executor.submit(_csr_matrix_to_anndata,
                                matrix_paths=cur_cov_paths,
                                obs_names=cur_ids,
                                chrom_size_path=chrom_size_path,
                                bin_size=bin_size,
                                mc_type=mc_type,
                                count_type=count_type,
                                step_size=step_size,
                                strandness=strandness,
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
                            mc_type=mc_type,
                            count_type=count_type,
                            step_size=step_size,
                            strandness=strandness,
                            output_path=output_prefix + f'{adata_chunk_id}.mc.h5ad',
                            compression=compression,
                            compression_opts=None)
            # cov
            executor.submit(_csr_matrix_to_anndata,
                            matrix_paths=cur_cov_paths,
                            obs_names=cur_ids,
                            chrom_size_path=chrom_size_path,
                            bin_size=bin_size,
                            mc_type=mc_type,
                            count_type=count_type,
                            step_size=step_size,
                            strandness=strandness,
                            output_path=output_prefix + f'{adata_chunk_id}.cov.h5ad',
                            compression=compression,
                            compression_opts=None)

    # remove temp file until everything finished
    subprocess.run(['rm', '-f'] + mc_path_list + cov_path_list)
    return


def _region_count_table_to_csr_npz(region_count_tables,
                                   region_id_map,
                                   output_prefix,
                                   compression=True):
    output_prefix = output_prefix.rstrip('.')
    mc_path = output_prefix + '.mc.npz'
    cov_path = output_prefix + '.cov.npz'
    if pathlib.Path(mc_path).exists() and pathlib.Path(cov_path).exists():
        return mc_path, cov_path

    count_table_path_list = parse_file_paths(region_count_tables)

    # read region_id to int_id (col location) map
    if isinstance(region_id_map, str):
        region_id_map = pd.read_msgpack(region_id_map).astype(np.uint32)
    else:
        region_id_map = region_id_map.astype(np.uint32)

    # determine matrix shape
    n_obj = len(count_table_path_list)
    n_feature = len(region_id_map)

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
                           dtype={'bin_id': str,
                                  'mc': np.uint32,
                                  'cov': np.uint32},
                           index_col=0)
        data.index = data.index.map(region_id_map)
        mc_matrix[obj_idx, data.index.values] = data['mc'].values
        cov_matrix[obj_idx, data.index.values] = data['cov'].values
    mc_matrix = mc_matrix.tocsr()
    ss.save_npz(mc_path, mc_matrix, compressed=compression)

    cov_matrix = cov_matrix.tocsr()
    ss.save_npz(cov_path, cov_matrix, compressed=compression)
    return mc_path, cov_path


def _csr_matrix_to_dataarray(matrix_table,
                             row_name, row_index,
                             col_name, col_index,
                             other_dim_info):
    total_mc_matrix = ss.vstack([ss.load_npz(path) for path in matrix_table['mc']]).todense()
    total_cov_matrix = ss.vstack([ss.load_npz(path) for path in matrix_table['cov']]).todense()

    data_arrays = []
    for count_type, matrix in zip(['mc', 'cov'], [total_mc_matrix, total_cov_matrix]):
        data_array = xr.DataArray(total_mc_matrix,
                                  dims=[row_name, col_name],
                                  coords={row_name: row_index,
                                          col_name: col_index.tolist()})
        _other_dim_info = other_dim_info.copy()
        _other_dim_info['count_type'] = count_type

        dim_axis = 2
        for dim_name, dim_coord in _other_dim_info.items():
            data_array = data_array.expand_dims(dim_name, axis=dim_axis)
            dim_axis += 1
            data_array.coords[dim_name] = [dim_coord]
        data_arrays.append(data_array)
    data_array = xr.concat(data_arrays, dim='count_type')
    return data_array


def aggregate_region_count_to_mcds(output_dir, dataset_name, chunk_size=100, row_name='cell'):
    # TODO write prepare mcds.py
    # TODO write test
    output_dir = pathlib.Path(output_dir)
    summary_df = pd.read_msgpack(output_dir / 'REGION_COUNT_SUMMARY.msg')
    file_uids = summary_df['file_id'].unique()

    csr_matrix_records = {}
    region_index_dict = {}
    for (mc_type, region_name, strandness), sub_summary_df in summary_df.groupby(
            ['mc_type', 'region_name', 'strandness']):
        if region_name not in region_index_dict:
            region_index = pd.read_msgpack(output_dir / f'REGION_ID_{region_name}.msg').index
            region_index.name = region_name
            region_index_dict[region_name] = region_index

        records = []
        for chunk_id, chunk_start in enumerate(range(0, file_uids.size, chunk_size)):
            file_id_chunk = file_uids[chunk_start: chunk_start + chunk_size]
            file_paths = sub_summary_df.loc[file_id_chunk]['file_path']

            mc_path, cov_path = _region_count_table_to_csr_npz(region_count_tables=file_paths,
                                                               region_id_map='./5E/REGION_ID_gene.msg',
                                                               output_prefix='5E/agg',
                                                               compression=True)
            records.append({'mc': mc_path, 'cov': cov_path})
        matrix_table = pd.DataFrame(records)
        if strandness.lower() == 'crick':
            strandness = '-'
        elif strandness.lower() == 'watson':
            strandness = '+'
        else:
            strandness = 'both'
        csr_matrix_records[(mc_type, region_name, strandness)] = matrix_table

    region_da_records = {}
    for region_name in summary_df['region_name'].unique():
        mc_type_da_records = []
        for mc_type in summary_df['mc_type'].unique():
            strand_da_records = []
            for strandness in summary_df['strandness'].unique():
                other_dim_info = {'mc_type': mc_type,
                                  'strand_type': strandness}
                matrix_table = csr_matrix_records[(mc_type, region_name, strandness)]

                dataarray = _csr_matrix_to_dataarray(matrix_table=matrix_table,
                                                     row_name=row_name,
                                                     row_index=file_uids,
                                                     col_name=region_name,
                                                     col_index=region_index_dict[region_name],
                                                     other_dim_info=other_dim_info)
                strand_da_records.append(dataarray)
            mc_type_da = xr.concat(strand_da_records)
            mc_type_da_records.append(mc_type_da)
        region_da_records[region_name] = xr.concat(mc_type_da_records)

    total_ds = xr.Dataset(region_da_records)
    total_ds.to_netcdf(output_dir / f'{dataset_name}.mcds')
    return
