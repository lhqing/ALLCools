import pathlib
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import msgpack
import pandas as pd
import xarray as xr

from .._doc import *
from .._open import open_allc
from ..utilities import parse_mc_pattern


def _count_data_to_xarray(total_data: dict, allc_series):
    """total_data: key is ALLC path, value is series of motif cov/mc count sum from all regions of that ALLC"""
    total_data = pd.DataFrame(total_data).reset_index()
    dir_data_list = []
    directions = []
    for direction, direction_df in total_data.groupby('level_0'):
        data_list = []
        mc_patterns = []
        for mc_pattern, pattern_df in direction_df.groupby('level_1'):
            pattern_df = pattern_df.iloc[:, 2:].set_index('level_2')
            pattern_df.index.name = 'motif'
            pattern_df.columns.name = 'cell'
            data_list.append(xr.DataArray(pattern_df))
            mc_patterns.append(mc_pattern)
        direction_da: xr.DataArray = xr.concat(data_list, dim='mc_type')
        direction_da.coords['mc_type'] = mc_patterns

        dir_data_list.append(direction_da)
        directions.append(direction)
    total_da: xr.DataArray = xr.concat(dir_data_list, dim='relative_strand')

    direction_transfer = {1: 'forward', 0: 'reverse'}
    total_da.coords['relative_strand'] = [direction_transfer[i] for i in directions]
    total_da.coords['cell'] = allc_series.index.tolist()
    return total_da


def _count_single_cmotif_bin_on_multiple_allc(cmotif_dict_path, allc_paths,
                                              region, count_binary, context_to_pattern):
    with open(cmotif_dict_path, 'rb') as f:
        c_motif_dict = msgpack.unpackb(f.read(), raw=True, use_list=False)

    mc_records = {}
    cov_records = {}
    for allc_path in allc_paths:
        with open_allc(allc_path, region=region) as allc:
            motif_mc_records = defaultdict(int)
            motif_cov_records = defaultdict(int)
            for line in allc:
                _, pos, _, context, mc, cov, _ = line.split('\t')
                pos = int(pos)
                try:
                    motifs = c_motif_dict[pos]
                    mc = int(mc)
                    cov = int(cov)
                    if count_binary:  # only for single cell, either methylated or unmethylated
                        for (rf_strand, motif_id) in motifs:
                            if mc > 0 and mc != cov:
                                # ambiguous mc and cov, just skip this
                                continue
                            motif_mc_records[
                                (rf_strand, context_to_pattern[context], motif_id)] += 0 if mc == 0 else 1
                            motif_cov_records[(rf_strand, context_to_pattern[context], motif_id)] += 1
                    else:
                        for (rf_strand, motif_id) in motifs:
                            motif_mc_records[(rf_strand, context_to_pattern[context], motif_id)] += mc
                            motif_cov_records[(rf_strand, context_to_pattern[context], motif_id)] += cov
                except KeyError:
                    continue
        mc_records[allc_path] = pd.Series(motif_mc_records)
        cov_records[allc_path] = pd.Series(motif_cov_records)
    return mc_records, cov_records


@doc_params(allc_table_doc=allc_table_doc,
            cpu_basic_doc=cpu_basic_doc)
def allc_motif_scan(allc_table, output_path, mc_contexts, c_motif_dir,
                    count_binary=True, cpu=1):
    """\
    Scan a list of ALLC files using a C-Motif database.
    C-Motif Database, can be generated via 'allcools generate-cmotif-database'
    Save the integrated multi-dimensional array into netCDF4 format using xarray.

    Parameters
    ----------
    allc_table
        {allc_table_doc}
    mc_contexts
    c_motif_dir
        A directory contains list of .msg files, each file records a dict,
        position is key, value is tuple of motif direction and id
    output_path
    count_binary
        Only use this for single cell allc, instead of sum mC or cov directly,
        will transfer mC and cov into [0, 1] when there is not ambiguity.
    cpu
        {cpu_basic_doc}

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

    c_motif_dir = pathlib.Path(c_motif_dir).absolute()
    allc_paths = allc_series.tolist()

    context_to_pattern = {}
    for pattern in mc_contexts:
        for context in parse_mc_pattern(pattern):
            if context in context_to_pattern:
                raise ValueError('Do not support overlapped pattern')
            else:
                context_to_pattern[context] = pattern

    lookup_table = pd.read_hdf(c_motif_dir / 'LOOKUP_TABLE.hdf')

    # TODO make this function more memory efficient

    region_mc_records = []
    region_cov_records = []
    with ProcessPoolExecutor(cpu) as executor:
        future_dict = {}
        for _, (chrom, start, end, file_name) in lookup_table.iterrows():
            region = f'{chrom}:{start}-{end}'
            cmotif_dict_path = str(c_motif_dir / file_name)
            future = executor.submit(_count_single_cmotif_bin_on_multiple_allc,
                                     cmotif_dict_path=cmotif_dict_path,
                                     allc_paths=allc_paths,
                                     region=region,
                                     count_binary=count_binary,
                                     context_to_pattern=context_to_pattern)
            future_dict[future] = region

        for future in as_completed(future_dict):
            region = future_dict[future]
            print(region)
            mc_records, cov_records = future.result()
            # the regions are not in order,
            # but this actually doesn't matter because we will sum them all together
            region_mc_records.append(mc_records)
            region_cov_records.append(cov_records)

    total_mc_data = {}
    total_cov_data = {}
    # aggregate regions for each ALLC
    for allc_path in allc_paths:
        # mc
        allc_mc_list = [mc_records[allc_path] for mc_records in region_mc_records]
        total_mc_data[allc_path] = pd.DataFrame(allc_mc_list).sum(axis=0)
        # cov
        allc_cov_list = [cov_records[allc_path] for cov_records in region_cov_records]
        total_cov_data[allc_path] = pd.DataFrame(allc_cov_list).sum(axis=0)

    # aggregate all ALLC
    total_mc_da = _count_data_to_xarray(total_mc_data, allc_series)
    total_cov_da = _count_data_to_xarray(total_cov_data, allc_series)
    total_da: xr.DataArray = xr.concat([total_mc_da, total_cov_da], dim='count_type')
    total_da.coords['count_type'] = ['mc', 'cov']

    with open(c_motif_dir / 'MOTIF_NAMES.msg', 'rb') as f:
        motif_names = msgpack.unpackb(f.read(), raw=False, use_list=False)
    motif_names = pd.Series({v: k for k, v in motif_names.items()})
    motif_names.index.name = 'motif'
    total_da.coords['motif_names'] = motif_names

    total_da = total_da.set_index(motif='motif_names')
    total_da.to_netcdf(output_path)
    return
