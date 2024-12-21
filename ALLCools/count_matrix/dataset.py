import pathlib
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from shutil import rmtree
import tempfile


import numpy as np
import pandas as pd
import pybedtools
import pysam
import xarray as xr
from numcodecs import blosc
from scipy import stats
import zarr, zarr.creation, zarr.convenience, zarr.hierarchy, zarr.storage

from ALLCools.utilities import parse_chrom_size, parse_mc_pattern

from .._doc import *

ALLOW_QUANT_TYPES = ["count", "hypo-score", "hyper-score"]


@lru_cache(99999)
def _bin_sf(cov, mc, p):
    if cov > mc:
        return stats.binom(cov, p).sf(mc)
    else:
        # cov == mc, sf = 0
        return 0


def _cell_sf(cell_count_df):
    mc_sum, cov_sum = cell_count_df.sum()
    p = mc_sum / (cov_sum + 0.000001)  # prevent empty allc error
    pv = cell_count_df.apply(lambda x: _bin_sf(x["cov"], x["mc"], p), axis=1).astype("float16")
    return pv


class _Quant:
    def __init__(self, mc_types, quant_type, kwargs):
        self.mc_types = mc_types
        self.quant_type = quant_type
        self.kwargs = kwargs


class _CountQuantifier:
    """Sum count by mC type for single region in single ALLC."""

    def __init__(self, mc_types):
        self.mc_types = mc_types
        self.mc_data = defaultdict(int)
        self.cov_data = defaultdict(int)
        return

    def read_line(self, line):
        """Read a line from an ALLC file."""
        chrom, pos, _, context, mc, cov, *_ = line.split("\t")
        self.mc_data[context] += int(mc)
        self.cov_data[context] += int(cov)
        return

    def summary(self):
        """Return a summary of the data."""
        mc_type_data = []
        for mc_type in self.mc_types:
            mc = sum(self.mc_data[k] for k in parse_mc_pattern(mc_type))
            cov = sum(self.cov_data[k] for k in parse_mc_pattern(mc_type))
            mc_type_data.append([mc, cov])
        return mc_type_data


def _determine_datasets(regions, quantifiers, chrom_size_path):
    """Determine datasets for each region."""
    tmpdir = tempfile.mkdtemp()
    chrom_sizes = parse_chrom_size(chrom_size_path)
    datasets = {}
    for pair in regions:
        if len(pair) != 2:
            raise ValueError(f"Can not understand {pair} in regions parameter: {regions}")
        name, region_path = pair
        # prepare regions
        if isinstance(region_path, (str, pathlib.Path)) and pathlib.Path(region_path).exists():
            region_bed_df = pd.read_csv(region_path, sep="\t", header=None)
            # remove additional chroms that do not occur in chrom_sizes
            region_bed_df = region_bed_df[region_bed_df.iloc[:, 0].isin(chrom_sizes)]
            # sort chroms
            region_bed_df = pybedtools.BedTool.from_dataframe(region_bed_df).sort(g=chrom_size_path).to_dataframe()

            if region_bed_df.shape[1] == 3:
                # add index
                print(
                    region_path,
                    "do not have index in its fourth column, adding it automatically. "
                    "If this is not desired, add a fourth column containing UNIQUE IDs to the BED file.",
                )
                region_bed_df[name] = (f"{name}_{i}" for i in range(region_bed_df.shape[0]))
            # check if name is unique()
            if region_bed_df.iloc[:, 3].duplicated().sum() > 0:
                raise ValueError(f"Region IDs in {region_path} (fourth column) are not unique.")
            # finally, set ID as index and only take the first three columns
            region_bed_df = region_bed_df.iloc[:, [0, 1, 2, 3]].set_index(region_bed_df.columns[3])
        else:
            try:
                # if region is int, generate chrom bins with bedtools
                region_size = int(region_path)
                region_bed_df = (
                    pybedtools.BedTool().makewindows(g=chrom_size_path, w=region_size, s=region_size).to_dataframe()
                )

                # set region index
                _dfs = []
                for chrom, chrom_df in region_bed_df.groupby("chrom"):
                    chrom_df = chrom_df.reset_index(drop=True)

                    def _id(i, c=chrom):
                        return f"{c}_{i}"

                    chrom_df.index = chrom_df.index.map(_id)
                    _dfs.append(chrom_df)
                region_bed_df = pd.concat(_dfs)

            except ValueError:
                raise ValueError(f"Can not understand region specification {region_path}")
        region_path = f"{tmpdir}/{name}.regions.csv"
        region_bed_df.to_csv(region_path)
        datasets[name] = {"regions": region_path, "quant": []}

    for quantifier in quantifiers:
        if len(quantifier) < 3:
            raise ValueError(
                f"Quantifier must have three parts, including "
                f'["NAME", "QUANT_TYPE", "MC_TYPE", "OTHER_KWARGS"], '
                f'where the "OTHER_KWARGS" are optional. '
                f"Got {quantifier}"
            )
        name, quant_type, mc_types, *other_kwargs = quantifier
        if name not in datasets:
            raise KeyError(f"Name {name} occur in quantifiers, but not found in regions.")
        kwargs = {}
        for kv in other_kwargs:
            k, v = kv.split("=")
            try:
                kwargs[k] = float(v)
            except ValueError:
                kwargs[k] = v

        # prepare mc_types
        mc_types = [i.strip() for i in mc_types.split(",")]
        # prepare quant_types

        if quant_type not in ALLOW_QUANT_TYPES:
            raise ValueError(f"QUANT_TYPE need to be in {ALLOW_QUANT_TYPES}, got {quant_type} in {quantifier}.")
        datasets[name]["quant"].append(_Quant(mc_types=mc_types, quant_type=quant_type, kwargs=kwargs))
    return datasets, tmpdir


def _count_single_region_set(allc_table, region_config, obs_dim, region_dim):
    """Get cell-by-region-by-mc_types count matrix, save to zarr."""
    total_mc_types = []
    for quant in region_config["quant"]:
        total_mc_types += quant.mc_types
    total_mc_types = list(set(total_mc_types))

    total_data = []
    for sample, allc_path in allc_table.items():
        with pysam.TabixFile(allc_path) as allc:
            region_ids = []
            sample_data = []
            region_chunks = pd.read_csv(region_config["regions"], index_col=0, chunksize=1000)
            for chunk in region_chunks:
                region_ids += chunk.index.tolist()
                for _, (chrom, start, end) in chunk.iterrows():
                    count_quant = _CountQuantifier(mc_types=total_mc_types)
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
                coords=[[sample], region_ids, total_mc_types, ["mc", "cov"]],
                dims=[obs_dim, region_dim, "mc_type", "count_type"]
            )
            total_data.append(data)
    total_data = xr.Dataset({f"{region_dim}_da": xr.concat(total_data, dim=obs_dim)})
    return total_data


def _calculate_pv(data, reverse_value, obs_dim, var_dim, cutoff=0.9):
    pv = []
    for cell in data.get_index(obs_dim):
        value = _cell_sf(data.sel(cell=cell).to_pandas())
        pv.append(value)
    pv = np.array(pv)

    if reverse_value:
        pv = 1 - pv

    # get rid of small values, save space and memory
    pv[pv < cutoff] = 0
    pv = xr.DataArray(pv, coords=[data.coords[obs_dim], data.coords[var_dim]], dims=[obs_dim, var_dim])
    pv = pv.astype("float16")
    return pv


def _count_single_zarr(
    allc_table, region_config, obs_dim, region_dim, chunk_start, regiongroup, count_dtype="uint32"
):
    """Process single region set and its quantifiers."""
    # count all ALLC and mC types that's needed for quantifiers if this region_dim
    count_ds = _count_single_region_set(
        allc_table=allc_table, region_config=region_config, obs_dim=obs_dim, region_dim=region_dim
    )

    # deal with count quantifiers
    count_mc_types = []
    for quant in region_config["quant"]:
        if quant.quant_type == "count":
            count_mc_types += quant.mc_types
    count_mc_types = list(set(count_mc_types))
    if len(count_mc_types) > 0:
        count_da = count_ds.sel(mc_type=count_mc_types)[f"{region_dim}_da"]
        max_int = np.iinfo(count_dtype).max
        count_da = xr.where(count_da > max_int, max_int, count_da)
        regiongroup[f"{region_dim}_da"][
            chunk_start : chunk_start + allc_table.index.size, :, :, :] = count_da.astype(count_dtype).data
    # deal with hypo-score, hyper-score quantifiers
    for quant in region_config["quant"]:
        if quant.quant_type == "hypo-score":
            for mc_type in quant.mc_types:
                data = _calculate_pv(
                    data=count_ds.sel(mc_type=mc_type)[f"{region_dim}_da"],
                    reverse_value=False,
                    obs_dim=obs_dim,
                    var_dim=region_dim,
                    **quant.kwargs,
                )
                regiongroup[f"{region_dim}_da_{mc_type}-hypo-score"][
                    chunk_start : chunk_start + allc_table.index.size, :
                ] = data.data
        elif quant.quant_type == "hyper-score":
            for mc_type in quant.mc_types:
                data = _calculate_pv(
                    count_ds.sel(mc_type=mc_type)[f"{region_dim}_da"],
                    reverse_value=True,
                    obs_dim=obs_dim,
                    var_dim=region_dim,
                    **quant.kwargs,
                )
                regiongroup[f"{region_dim}_da_{mc_type}-hyper-score"][chunk_start : chunk_start + allc_table.index.size, :] = data.data

    return True


@doc_params(
    generate_dataset_doc=generate_dataset_doc,
    allc_table_doc=allc_table_doc,
    chrom_size_path_doc=chrom_size_path_doc,
    regions_doc=generate_dataset_regions_doc,
    quantifiers_doc=generate_dataset_quantifiers_doc,
    obs_dim_doc=generate_dataset_obs_dim_doc,
    cpu_basic_doc=cpu_basic_doc,
    chunk_size_doc=generate_dataset_chunk_size_doc,
)
def generate_dataset(
    allc_table, output_path, regions, quantifiers, chrom_size_path, obs_dim="cell", cpu=1, chunk_size=None
):
    """\
    {generate_dataset_doc}

    Parameters
    ----------
    allc_table
        {allc_table_doc}
    output_path
        Output path of the MCDS dataset
    regions
        {regions_doc}
    quantifiers
        {quantifiers_doc}
    chrom_size_path
        {chrom_size_path_doc}
    obs_dim
        {obs_dim_doc}
    cpu
        {cpu_basic_doc}
    chunk_size
        {chunk_size_doc}

    Returns
    -------
    output_path
    """
    if isinstance(allc_table, (str, pathlib.Path)):
        allc_table = pd.read_csv(allc_table, sep="\t", header=None, index_col=0).squeeze()
        allc_table.index.name = obs_dim

    # determine index length and str dtype
    max_length = allc_table.index.map(lambda idx: len(idx)).max()

    # determine parallel chunk size
    n_sample = allc_table.size
    if chunk_size is None:
        chunk_size = min(n_sample, 50)

    # prepare regions and determine quantifiers
    pathlib.Path(output_path).mkdir(exist_ok=True)
    z = zarr.storage.DirectoryStore(path=output_path)
    root = zarr.hierarchy.group(store = z, overwrite = True)
    datasets, tmpdir = _determine_datasets(regions, quantifiers, chrom_size_path)
    # copy chrom_size_path to output_path
    subprocess.run(["cp", "-f", chrom_size_path, f"{output_path}/chrom_sizes.txt"], check=True)
    for region_dim, region_config in datasets.items():
        regiongroup = root.create_group(region_dim)
        # save region coords to the ds
        bed = pd.read_csv(f"{tmpdir}/{region_dim}.regions.csv", index_col=0)
        bed.columns = [f"{region_dim}_chrom", f"{region_dim}_start", f"{region_dim}_end"]
        bed.index.name = region_dim
        region_size = bed.index.size
        dsobs = regiongroup.array(
            name=obs_dim,
            data=allc_table.index.values,
            chunks=(chunk_size),
            dtype=f"<U{max_length}"
        )
        dsobs.attrs['_ARRAY_DIMENSIONS'] = [obs_dim]
        # append region bed to the saved ds
        ds = xr.Dataset()
        for col, data in bed.items():
            ds.coords[col] = data
        ds.coords[region_dim] = bed.index.values
        # change object dtype to string
        for k in ds.coords.keys():
            if ds.coords[k].dtype == "O":
                ds.coords[k] = ds.coords[k].astype(str)
        ds.to_zarr(f"{output_path}/{region_dim}", mode="w")
        count_mc_types = []
        for quant in region_config["quant"]:
            if quant.quant_type == "count":
                count_mc_types += quant.mc_types
        count_mc_types = list(set(count_mc_types))
        if len(count_mc_types) > 0:
            DA = regiongroup.empty(
                name=f"{region_dim}_da",
                shape=(n_sample, region_size, len(count_mc_types), 2),
                chunks=(chunk_size, region_size, len(count_mc_types), 2),
                dtype="uint32"
            )
            DA.attrs['_ARRAY_DIMENSIONS']=[obs_dim, region_dim, "mc_type", "count_type"]
            count = regiongroup.array(
                name="count_type",
                data=(["mc", "cov"]),
                dtype="<U3"
            )
            count.attrs['_ARRAY_DIMENSIONS']=["count_type"]
            mc = regiongroup.array(
                name="mc_type",
                data=count_mc_types,
                dtype="<U3"
            )
            mc.attrs['_ARRAY_DIMENSIONS']=["mc_type"]
        # deal with hypo-score, hyper-score quantifiers
        for quant in region_config["quant"]:
            if quant.quant_type == "hypo-score":
                for mc_type in quant.mc_types:
                    hypo = regiongroup.empty (
                        name = f"{region_dim}_da_{mc_type}-hypo-score",
                        shape=(allc_table.size, region_size),
                        chunks = (chunk_size, region_size),
                        dtype = "float16"
                    )
                    hypo.attrs['_ARRAY_DIMENSIONS']=[obs_dim, region_dim]
            elif quant.quant_type == "hyper-score":
                for mc_type in quant.mc_types:
                    hyper = regiongroup.empty (
                        name = f"{region_dim}_da_{mc_type}-hyper-score",
                        shape=(allc_table.size, region_size),
                        chunks = (chunk_size, region_size),
                        dtype = "float16"
                    )
                    hyper.attrs['_ARRAY_DIMENSIONS']=[obs_dim, region_dim]
    blosc.use_threads = False
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        # parallel on allc chunks and region_sets levels
        for i, chunk_start in enumerate(range(0, n_sample, chunk_size)):
            allc_chunk = allc_table[chunk_start : chunk_start + chunk_size]
            for region_dim, region_config in datasets.items():
                f = exe.submit(
                    _count_single_zarr,
                    allc_table=allc_chunk,
                    region_config=region_config,
                    obs_dim=obs_dim,
                    region_dim=region_dim,
                    chunk_start=chunk_start,
                    regiongroup=regiongroup,
                )
                futures[f] = (region_dim, i)

        for f in as_completed(futures):
            region_dim, i = futures[f]
            print(f"Chunk {i} of {region_dim} returned")
    blosc.use_threads = None
    from ..mcds.utilities import update_dataset_config

    update_dataset_config(
        output_path,
        config={
            "region_dim": None,
            "ds_region_dim": {region_dim: region_dim for region_dim in datasets.keys()},
            "ds_sample_dim": {region_dim: obs_dim for region_dim in datasets.keys()},
        },
    )
    zarr.convenience.consolidate_metadata(z)
    return output_path
