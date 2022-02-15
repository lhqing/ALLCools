from ALLCools.utilities import parse_chrom_size, parse_mc_pattern
import pybedtools
from collections import defaultdict
import pathlib
import pandas as pd
import xarray as xr
import numpy as np
import pysam
import subprocess
from functools import lru_cache
from scipy import stats
from concurrent.futures import ProcessPoolExecutor, as_completed
from shutil import rmtree

ALLOW_QUANT_TYPES = ['count', 'hypo-score', 'hyper-score']


@lru_cache(99999)
def bin_sf(cov, mc, p):
    if cov > mc:
        return stats.binom(cov, p).sf(mc)
    else:
        # cov == mc, sf = 0
        return 0


def cell_sf(cell_count_df):
    mc_sum, cov_sum = cell_count_df.sum()
    p = mc_sum / (cov_sum + 0.000001)  # prevent empty allc error
    pv = cell_count_df.apply(lambda x: bin_sf(x["cov"], x["mc"], p), axis=1).astype(
        "float16"
    )
    return pv


class Quant:
    def __init__(self, mc_types, quant_type, kwargs):
        self.mc_types = mc_types
        self.quant_type = quant_type
        self.kwargs = kwargs


class CountQuantifier:
    """Sum count by mC type for single region in single ALLC"""

    def __init__(self, mc_types):
        self.mc_types = mc_types
        self.mc_data = defaultdict(int)
        self.cov_data = defaultdict(int)
        return

    def read_line(self, line):
        chrom, pos, _, context, mc, cov, *_ = line.split("\t")
        self.mc_data[context] += int(mc)
        self.cov_data[context] += int(cov)
        return

    def summary(self):
        mc_type_data = []
        for mc_type in self.mc_types:
            mc = sum((self.mc_data[k] for k in parse_mc_pattern(mc_type)))
            cov = sum((self.cov_data[k] for k in parse_mc_pattern(mc_type)))
            mc_type_data.append([mc, cov])
        return mc_type_data


def determine_datasets(regions, quantifiers, chrom_size_path, tmp_dir):
    tmp_dir = pathlib.Path(tmp_dir).absolute()
    tmp_dir.mkdir(exist_ok=True, parents=True)

    chrom_sizes = parse_chrom_size(chrom_size_path)
    datasets = {}
    for pair in regions:
        if len(pair) != 2:
            raise ValueError(f'Can not understand {pair} in regions parameter: {regions}')
        name, region_path = pair
        # prepare regions
        if isinstance(region_path, (str, pathlib.Path)) \
                and pathlib.Path(region_path).exists():
            region_bed_df = pd.read_csv(region_path, sep='\t', header=None)
            # remove additional chroms that do not occur in chrom_sizes
            region_bed_df = region_bed_df[region_bed_df.iloc[:, 0].isin(
                chrom_sizes)]
            # sort chroms
            region_bed_df = pybedtools.BedTool.from_dataframe(
                region_bed_df).sort(g=chrom_size_path).to_dataframe()

            if region_bed_df.shape[1] == 3:
                # add index
                print(
                    region_path,
                    'do not have index in its fourth column, adding it automatically. '
                    'If this is not desired, add a fourth column containing UNIQUE IDs to the BED file.'
                )
                region_bed_df[name] = (f'{name}_{i}'
                                       for i in range(region_bed_df.shape[0]))
            # check if name is unique()
            if region_bed_df.iloc[:, 3].duplicated().sum() > 0:
                raise ValueError(
                    f'Region IDs in {region_path} (fourth column) are not unique.'
                )
            # finally, set ID as index and only take the first three columns
            region_bed_df = region_bed_df.iloc[:, [0, 1, 2, 3]].set_index(
                region_bed_df.columns[3])
        else:
            try:
                # if region is int, generate chrom bins with bedtools
                region_size = int(region_path)
                region_bed_df = pybedtools.BedTool().makewindows(
                    g=chrom_size_path, w=region_size,
                    s=region_size).to_dataframe()

                # set region index
                _dfs = []
                for chrom, chrom_df in region_bed_df.groupby('chrom'):
                    chrom_df = chrom_df.reset_index(drop=True)
                    chrom_df.index = chrom_df.index.map(lambda i: f'{chrom}_{i}')
                    _dfs.append(chrom_df)
                region_bed_df = pd.concat(_dfs)

            except ValueError:
                raise ValueError(
                    f'Can not understand region specification {region_path}')
        region_path = f'{tmp_dir}/{name}.regions.csv'
        region_bed_df.to_csv(region_path)
        datasets[name] = {'regions': region_path, 'quant': []}

    for quantifier in quantifiers:
        if len(quantifier) < 3:
            raise ValueError(
                f'Quantifier must have three parts, including '
                f'["NAME", "QUANT_TYPE", "MC_TYPE", "OTHER_KWARGS"], '
                f'where the "OTHER_KWARGS" are optional. '
                f'Got {quantifier}')
        name, quant_type, mc_types, *other_kwargs = quantifier
        if name not in datasets:
            raise KeyError(
                f'Name {name} occur in quantifiers, but not found in regions.')
        kwargs = {}
        for kv in other_kwargs:
            k, v = kv.split('=')
            try:
                kwargs[k] = float(v)
            except ValueError:
                kwargs[k] = v

        # prepare mc_types
        mc_types = [i.strip() for i in mc_types.split(',')]
        # prepare quant_types

        if quant_type not in ALLOW_QUANT_TYPES:
            raise ValueError(
                f'QUANT_TYPE need to be in {ALLOW_QUANT_TYPES}, got {quant_type} in {quantifier}.'
            )
        datasets[name]['quant'].append(
            Quant(mc_types=mc_types, quant_type=quant_type, kwargs=kwargs))
    return datasets


def count_single_region_set(allc_table, region_config, obs_dim, region_dim):
    """Get cell-by-region-by-mc_types count matrix, save to zarr"""
    total_mc_types = []
    for quant in region_config['quant']:
        total_mc_types += quant.mc_types
    total_mc_types = list(set(total_mc_types))

    total_data = []
    for sample, allc_path in allc_table.items():
        with pysam.TabixFile(allc_path) as allc:
            region_ids = []
            sample_data = []
            region_chunks = pd.read_csv(region_config['regions'],
                                        index_col=0,
                                        chunksize=1000)
            for chunk in region_chunks:
                region_ids += chunk.index.tolist()
                for region, (chrom, start, end) in chunk.iterrows():
                    count_quant = CountQuantifier(mc_types=total_mc_types)
                    try:
                        allc_lines = allc.fetch(chrom, start, end)
                        for line in allc_lines:
                            count_quant.read_line(line)
                    except ValueError:
                        # got value error, this chrom not exist in allc
                        pass
                    sample_data.append(count_quant.summary())
            data = xr.DataArray(
                np.array([sample_data]),
                coords=[[sample], region_ids, total_mc_types, ['mc', 'cov']],
                dims=[obs_dim, region_dim, 'mc_type', 'count_type'])
            total_data.append(data)
    total_data = xr.Dataset(
        {f'{region_dim}_da': xr.concat(total_data, dim=obs_dim)})
    return total_data


def calculate_pv(data, reverse_value, obs_dim, var_dim, cutoff=0.9):
    pv = []
    for cell in data.get_index(obs_dim):
        value = cell_sf(data.sel(cell=cell).to_pandas())
        pv.append(value)
    pv = np.array(pv)

    if reverse_value:
        pv = 1 - pv

    # get rid of small values, save space and memory
    pv[pv < cutoff] = 0
    pv = xr.DataArray(pv,
                      coords=[data.coords[obs_dim], data.coords[var_dim]],
                      dims=[obs_dim, var_dim])
    pv = pv.astype('float16')
    return pv


def count_single_zarr(allc_table,
                      region_config,
                      obs_dim,
                      region_dim,
                      output_path,
                      obs_dim_dtype,
                      count_dtype='uint32'):
    """process single region set and its quantifiers"""
    # count all ALLC and mC types that's needed for quantifiers if this region_dim
    count_ds = count_single_region_set(allc_table=allc_table,
                                       region_config=region_config,
                                       obs_dim=obs_dim,
                                       region_dim=region_dim)

    total_ds = {}
    # deal with count quantifiers
    count_mc_types = []
    for quant in region_config['quant']:
        if quant.quant_type == 'count':
            count_mc_types += quant.mc_types
    count_mc_types = list(set(count_mc_types))
    if len(count_mc_types) > 0:
        count_da = count_ds.sel(mc_type=count_mc_types)[f'{region_dim}_da']
        max_int = np.iinfo(count_dtype).max
        count_da = xr.where(count_da > max_int, max_int, count_da)
        total_ds[f'{region_dim}_da'] = count_da.astype(count_dtype)

    # deal with hypo-score, hyper-score quantifiers
    for quant in region_config['quant']:
        if quant.quant_type == 'hypo-score':
            for mc_type in quant.mc_types:
                data = calculate_pv(
                    data=count_ds.sel(mc_type=mc_type)[f'{region_dim}_da'],
                    reverse_value=False,
                    obs_dim=obs_dim,
                    var_dim=region_dim,
                    **quant.kwargs)
                total_ds[f'{region_dim}_da_{mc_type}-hypo-score'] = data
        elif quant.quant_type == 'hyper-score':
            for mc_type in quant.mc_types:
                data = calculate_pv(
                    count_ds.sel(mc_type=mc_type)[f'{region_dim}_da'],
                    reverse_value=True,
                    obs_dim=obs_dim,
                    var_dim=region_dim,
                    **quant.kwargs)
                total_ds[f'{region_dim}_da_{mc_type}-hyper-score'] = data
    total_ds = xr.Dataset(total_ds)
    total_ds.coords[obs_dim] = total_ds.coords[obs_dim].astype(obs_dim_dtype)
    total_ds.to_zarr(output_path, mode='w')
    return output_path


def generate_dataset(allc_table, output_path, regions, quantifiers, chrom_size_path,
                     obs_dim='cell', cpu=1, chunk_size=None):
    """
    Generate multiple methylation datasets with a set of allc_table,
    a list of region sets and quantifiers for each region set.

    Parameters
    ----------
    allc_table
    output_path
    regions
    quantifiers
    chrom_size_path
    obs_dim
    cpu
    chunk_size

    Returns
    -------
    output_path
    """
    if isinstance(allc_table, (str, pathlib.Path)):
        allc_table = pd.read_csv(allc_table,
                                 sep='\t',
                                 header=None,
                                 index_col=0,
                                 squeeze=True)
        allc_table.index.name = obs_dim

    # determine index length and str dtype
    max_length = allc_table.index.map(lambda idx: len(idx)).max()
    obs_dim_dtype = f'<U{max_length}'

    # determine parallel chunk size
    n_sample = allc_table.size
    if chunk_size is None:
        chunk_size = min(n_sample, 50)

    # prepare regions and determine quantifiers
    pathlib.Path(output_path).mkdir(exist_ok=True)
    tmp_dir = f'{output_path}_tmp'
    datasets = determine_datasets(regions, quantifiers, chrom_size_path, tmp_dir)

    # copy chrom_size_path to output_path
    subprocess.run(['cp', '-f', chrom_size_path,
                    f'{output_path}/chrom_sizes.txt'],
                   check=True)

    chunk_records = defaultdict(dict)
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        # parallel on allc chunks and region_sets levels
        for i, chunk_start in enumerate(range(0, n_sample, chunk_size)):
            allc_chunk = allc_table[chunk_start:chunk_start + chunk_size]
            for region_dim, region_config in datasets.items():
                chunk_path = f'{tmp_dir}/chunk_{region_dim}_{chunk_start}.zarr'
                f = exe.submit(count_single_zarr,
                               allc_table=allc_chunk,
                               region_config=region_config,
                               obs_dim=obs_dim,
                               region_dim=region_dim,
                               output_path=chunk_path,
                               obs_dim_dtype=obs_dim_dtype)
                futures[f] = (region_dim, i)

        for f in as_completed(futures):
            region_dim, i = futures[f]
            chunk_path = f.result()
            print(f'Chunk {i} of {region_dim} returned')
            chunk_records[region_dim][i] = chunk_path

    for region_dim, chunks in chunk_records.items():
        # write chunk in order
        chunk_paths = pd.Series(chunks).sort_index().tolist()
        for i, chunk_path in enumerate(chunk_paths):
            ds = xr.open_zarr(chunk_path).load()
            # dump chunk to final place
            if i == 0:
                # first chunk
                ds.to_zarr(f'{output_path}/{region_dim}',
                           mode='w')
            else:
                # append
                ds.to_zarr(f'{output_path}/{region_dim}',
                           append_dim=obs_dim)
            rmtree(chunk_path)

        # save region coords to the ds
        bed = pd.read_csv(f'{tmp_dir}/{region_dim}.regions.csv', index_col=0)
        bed.columns = [f'{region_dim}_chrom', f'{region_dim}_start', f'{region_dim}_end']
        bed.index.name = region_dim
        # append region bed to the saved ds
        ds = xr.Dataset()
        for col, data in bed.items():
            ds.coords[col] = data
        # change object dtype to string
        for k in ds.coords.keys():
            if ds.coords[k].dtype == 'O':
                ds.coords[k] = ds.coords[k].astype(str)
        ds.to_zarr(f'{output_path}/{region_dim}', mode='a')

    # delete tmp
    rmtree(tmp_dir)

    from ..mcds.utilities import update_dataset_config
    update_dataset_config(output_path,
                          config={"region_dim": None,
                                  "ds_region_dim": {region_dim: region_dim
                                                    for region_dim in datasets.keys()},
                                  "ds_sample_dim": {region_dim: obs_dim
                                                    for region_dim in datasets.keys()}})
    return output_path
