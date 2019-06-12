from collections import defaultdict

import msgpack
import pandas as pd
import xarray as xr
import pathlib

from .._open import open_allc
from ..utilities import parse_mc_pattern, parse_file_paths


def count_data_to_xarray(total_data: dict):
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
    total_da: xr.DataArray = xr.concat(dir_data_list, dim='direction')
    total_da.coords['direction'] = directions
    return total_da


def allc_count_motif(allc_paths, output_path, mc_contexts, c_motif_dir, count_binary=True):
    """

    Parameters
    ----------
    allc_paths
    mc_contexts
    c_motif_dir
        A directory contains list of .msg files, each file records a dict,
        position is key, value is tuple of motif direction and id
    output_path
    count_binary
        Only use this for single cell allc

    Returns
    -------

    """
    c_motif_dir = pathlib.Path(c_motif_dir).absolute()
    allc_paths = parse_file_paths(allc_paths)

    context_to_pattern = {}
    for pattern in mc_contexts:
        for context in parse_mc_pattern(pattern):
            if context in context_to_pattern:
                raise ValueError('Do not support overlapped pattern')
            else:
                context_to_pattern[context] = pattern

    lookup_table = pd.read_msgpack(c_motif_dir / 'lookup_table.msg')
    allc_mc_records = defaultdict(list)
    allc_cov_records = defaultdict(list)
    for _, (chrom, start, end, file_name) in lookup_table.iterrows():
        region = f'{chrom}:{start}-{end}'
        print(region)
        with open(c_motif_dir / file_name, 'rb') as f:
            c_motif_dict = msgpack.unpackb(f.read(), raw=True, use_list=False)
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
                        if count_binary:  # only for single cell
                            for (rf_strand, motif_id) in motifs:
                                motif_mc_records[
                                    (rf_strand, context_to_pattern[context], motif_id)] += 0 if mc == 0 else 1
                                motif_cov_records[(rf_strand, context_to_pattern[context], motif_id)] += 1
                        else:
                            for (rf_strand, motif_id) in motifs:
                                motif_mc_records[(rf_strand, context_to_pattern[context], motif_id)] += mc
                                motif_cov_records[(rf_strand, context_to_pattern[context], motif_id)] += cov
                    except KeyError:
                        continue
            allc_mc_records[allc_path].append(pd.Series(motif_mc_records))
            allc_cov_records[allc_path].append(pd.Series(motif_cov_records))

    total_mc_data = {allc_path: pd.DataFrame(data).sum(axis=0) for allc_path, data in allc_mc_records.items()}
    total_cov_data = {allc_path: pd.DataFrame(data).sum(axis=0) for allc_path, data in allc_cov_records.items()}
    total_mc_da = count_data_to_xarray(total_mc_data)
    total_cov_da = count_data_to_xarray(total_cov_data)
    total_da: xr.DataArray = xr.concat([total_mc_da, total_cov_da], dim='count_type')
    total_da.coords['count_type'] = ['mc', 'cov']

    with open(c_motif_dir / 'motif_names.msg', 'rb') as f:
        motif_names = msgpack.unpackb(f.read(), raw=False, use_list=False)
    motif_names = pd.Series({v: k for k, v in motif_names.items()})
    motif_names.index.name = 'motif'
    total_da.coords['motif_names'] = motif_names

    total_da.to_netcdf(output_path)
    return
