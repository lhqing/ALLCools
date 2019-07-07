import pathlib
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from math import ceil

import numpy as np
import pandas as pd
import scipy.sparse as ss
import xarray as xr

from .._allc_to_region_count import batch_allc_to_region_count
from .._doc import *
from ..schema.mcds_schema import *
from ..utilities import parse_file_paths


def _region_count_table_to_csr_npz(region_count_tables,
                                   region_id_map,
                                   output_prefix,
                                   compression=True):
    """
    helper func of _aggregate_region_count_to_mcds

    Take a list of region count table paths, read,
    aggregate them into a 2D sparse matrix and save the mC and COV separately.
    This function don't take care of any path selection, but assume all region_count_table is homogeneous type
    It return the saved file path
    """
    output_prefix = output_prefix.rstrip('.')
    mc_path = output_prefix + '.mc.npz'
    cov_path = output_prefix + '.cov.npz'
    if pathlib.Path(mc_path).exists() and pathlib.Path(cov_path).exists():
        return mc_path, cov_path

    if len(region_count_tables) > 1:
        count_table_path_list = parse_file_paths(region_count_tables)
    else:
        # only one path in the list, which is possible for last chunk
        count_table_path_list = [pathlib.Path(region_count_tables[0]).resolve()]

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
        # TODO save internal region count table as msg
        data = pd.read_csv(file_path,
                           usecols=[3, 4, 5],
                           sep='\t',
                           header=None,
                           names=['bin_id', 'mc', 'cov'],
                           dtype={'bin_id': str,
                                  'mc': np.uint32,
                                  'cov': np.uint32},
                           index_col=0)
        if data.shape[0] == 0:
            # for rare case where the site_bed is empty
            continue
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
    """
    helper func of _aggregate_region_count_to_mcds

    This function aggregate sparse array files into a single xarray.DataArray,
    combining cell chunks, mc/cov count type together.
    The matrix_table provide all file paths, each row is for a cell chunk, with mc and cov matrix path separately.
    """
    total_mc_matrix = ss.vstack([ss.load_npz(path) for path in matrix_table['mc']]).todense()
    total_cov_matrix = ss.vstack([ss.load_npz(path) for path in matrix_table['cov']]).todense()

    data_arrays = []
    for count_type, matrix in zip(['mc', 'cov'], [total_mc_matrix, total_cov_matrix]):
        data_array = xr.DataArray(matrix,
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


def _aggregate_region_count_to_mcds(output_dir,
                                    dataset_name,
                                    chunk_size=100,
                                    row_name='cell',
                                    cpu=1):
    """
    This function aggregate all the region count table into a single mcds
    """
    # TODO write test
    output_dir = pathlib.Path(output_dir)
    summary_df = pd.read_msgpack(output_dir / 'REGION_COUNT_SUMMARY.msg')
    file_uids = summary_df['file_id'].unique()

    region_index_dict = {}
    additional_coords = {}
    with ProcessPoolExecutor(cpu) as executor:
        # aggregate count table into 2D sparse array, parallel in cell chunks
        future_dict = {}
        for (mc_type, region_name, strandness), sub_summary_df in summary_df.groupby(
                ['mc_type', 'region_name', 'strandness']):
            sub_summary_df = sub_summary_df.set_index('file_id')
            if region_name not in region_index_dict:
                region_index = pd.read_msgpack(output_dir / f'REGION_ID_{region_name}.msg').index
                region_index.name = region_name
                region_index_dict[region_name] = region_index

                region_bed = pd.read_msgpack(output_dir / f'REGION_BED_{region_name}.msg')
                for col, value in region_bed.iteritems():
                    _col = f'{region_name}_{col}'
                    value.index.name = region_name
                    additional_coords[_col] = value

            for chunk_id, chunk_start in enumerate(range(0, file_uids.size, chunk_size)):
                file_id_chunk = file_uids[chunk_start: chunk_start + chunk_size]
                file_paths = sub_summary_df.loc[file_id_chunk]['file_path'].tolist()
                future = executor.submit(_region_count_table_to_csr_npz,
                                         region_count_tables=file_paths,
                                         region_id_map=str(output_dir / f'REGION_ID_{region_name}.msg'),
                                         output_prefix=str(
                                             output_dir / f'{dataset_name}_{region_name}_'
                                             f'{mc_type}_{strandness}_{chunk_id}'),
                                         compression=True)
                future_dict[future] = (mc_type, region_name, strandness, chunk_id)

        records = defaultdict(list)
        for future in as_completed(future_dict):
            mc_type, region_name, strandness, chunk_id = future_dict[future]
            try:
                mc_path, cov_path = future.result()
                records[(mc_type, region_name, strandness)] \
                    .append({'mc': mc_path, 'cov': cov_path, 'index': chunk_id})
            except Exception as e:
                print(f'Error when calculating mc-{mc_type} region-{region_name} '
                      f'strand-{strandness} chunk-{chunk_id}.')
                raise e

        # IMPORTANT order csr_matrix_records by chunk_id
        csr_matrix_records = {k: pd.DataFrame(v).set_index('index').sort_index()
                              for k, v in records.items()}

    # aggregate all the sparse array into a single mcds, this step load everything and need large memory
    region_da_records = {}
    for region_name in summary_df['region_name'].unique():
        mc_type_da_records = []
        for mc_type in summary_df['mc_type'].unique():
            strand_da_records = []
            for strandness in summary_df['strandness'].unique():
                matrix_table = csr_matrix_records[(mc_type, region_name, strandness)]
                if strandness.lower() == 'crick':
                    strandness = '-'
                elif strandness.lower() == 'watson':
                    strandness = '+'
                else:
                    strandness = 'both'
                other_dim_info = {'mc_type': mc_type,
                                  'strand_type': strandness}

                dataarray = _csr_matrix_to_dataarray(matrix_table=matrix_table,
                                                     row_name=row_name,
                                                     row_index=file_uids,
                                                     col_name=region_name,
                                                     col_index=region_index_dict[region_name],
                                                     other_dim_info=other_dim_info)
                strand_da_records.append(dataarray)
            mc_type_da = xr.concat(strand_da_records, dim='strand_type')
            mc_type_da_records.append(mc_type_da)
        region_da_records[region_name + "_da"] = xr.concat(mc_type_da_records, dim='mc_type')

    total_ds = xr.Dataset(region_da_records)
    total_ds.coords.update(additional_coords)

    return total_ds


@doc_params(chrom_size_path_doc=chrom_size_path_doc,
            mc_contexts_doc=mc_contexts_doc,
            cov_cutoff_doc=cov_cutoff_doc,
            split_strand_doc=split_strand_doc,
            bin_sizes_doc=bin_sizes_doc,
            cpu_basic_doc=cpu_basic_doc,
            region_bed_paths_doc=region_bed_paths_doc,
            region_bed_names_doc=region_bed_names_doc)
def generate_mcds(allc_table,
                  output_prefix,
                  chrom_size_path,
                  mc_contexts,
                  split_strand=False,
                  bin_sizes=None,
                  region_bed_paths=None,
                  region_bed_names=None,
                  cov_cutoff=9999,
                  cpu=1,
                  remove_tmp=True,
                  max_per_mcds=3072,
                  cell_chunk_size=100):
    """\
    Generate MCDS from a list of ALLC file provided with file id.

    Parameters
    ----------
    allc_table
        Contain all the ALLC file information in 2 columns: 1. file_uid, 2. file_path. No header
    output_prefix
        Output prefix of the MCDS
    chrom_size_path
        {chrom_size_path_doc}
    mc_contexts
        {mc_contexts_doc}
    split_strand
        {split_strand_doc}
    bin_sizes
        {bin_sizes_doc}
    region_bed_paths
        {region_bed_paths_doc}
    region_bed_names
        {region_bed_names_doc}
    cov_cutoff
        {cov_cutoff_doc}
    cpu
        {cpu_basic_doc}
    remove_tmp
        Whether to remove the temp directory for generating MCDS
    max_per_mcds
        Maximum number of ALLC files to aggregate into 1 MCDS, if number of ALLC provided > max_per_mcds,
        will generate MCDS in chunks, with same prefix provided.
    cell_chunk_size
        Size of cell chunk in parallel aggregation. Do not have any effect on results.
        Large chunksize needs large memory.

    Returns
    -------

    """
    if isinstance(allc_table, str):
        allc_series = pd.read_csv(allc_table, header=None, index_col=0, squeeze=True, sep='\t')
        if not isinstance(allc_series, pd.Series):
            raise ValueError('allc_table malformed, should only have 2 columns, 1. file_uid, 2. file_path')
    else:
        allc_series = allc_table
    if allc_series.index.duplicated().sum() != 0:
        raise ValueError('allc_table file uid have duplicates (1st column)')

    # if allc files exceed max_per_mcds, save them into chunks
    if allc_series.size > max_per_mcds:
        mcds_n_chunk = ceil(allc_series.size / max_per_mcds)
        chunk_size = ceil(allc_series.size / mcds_n_chunk)
        allc_series_chunks = [allc_series[chunk_start:chunk_start + chunk_size]
                              for chunk_start in range(0, allc_series.size, chunk_size)]
        print(f'Number of ALLC files {allc_series.size} > max_per_mcds {max_per_mcds}, ')
        print(f'will generate MCDS in {len(allc_series_chunks)} chunks.')

        # mcds chunks execute sequentially
        for chunk_id, allc_series_chunk in enumerate(allc_series_chunks):
            generate_mcds(allc_table=allc_series_chunk,
                          output_prefix=f'{output_prefix}_{chunk_id}',
                          chrom_size_path=chrom_size_path,
                          mc_contexts=mc_contexts,
                          split_strand=split_strand,
                          bin_sizes=bin_sizes,
                          region_bed_paths=region_bed_paths,
                          region_bed_names=region_bed_names,
                          cov_cutoff=cov_cutoff,
                          cpu=cpu,
                          remove_tmp=remove_tmp,
                          max_per_mcds=max_per_mcds,
                          cell_chunk_size=cell_chunk_size)
        return

    # check region bed names
    if region_bed_names is not None:
        for region_bed_name in region_bed_names:
            check_custom_dim_name(region_bed_name)

    output_prefix = output_prefix.rstrip('.')
    output_dir = pathlib.Path(output_prefix + '.tmp_dir')
    output_dir.mkdir()
    dataset_name = pathlib.Path(output_prefix).name

    # count all the ALLC files, generate region count tables in a temp dir
    print('Count ALLC files')
    batch_allc_to_region_count(allc_series=allc_series,
                               output_dir=output_dir,
                               chrom_size_path=chrom_size_path,
                               mc_contexts=mc_contexts,
                               split_strand=split_strand,
                               bin_sizes=bin_sizes,
                               region_bed_paths=region_bed_paths,
                               region_bed_names=region_bed_names,
                               cov_cutoff=cov_cutoff,
                               cpu=cpu)

    # aggregate all the region count table into a mcds
    print('Aggregate Region Count Table')
    total_ds = _aggregate_region_count_to_mcds(output_dir=output_dir,
                                               dataset_name=dataset_name,
                                               chunk_size=cell_chunk_size,
                                               row_name='cell',
                                               cpu=cpu)
    if not output_prefix.endswith('.mcds'):
        output_path = output_prefix + '.mcds'
    else:
        output_path = output_prefix
    total_ds.to_netcdf(output_path)

    # remove temp dir
    if remove_tmp:
        subprocess.run(['rm', '-rf', str(output_dir)], check=True)
    return
