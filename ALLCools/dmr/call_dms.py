import numpy as np
import pathlib
import pysam
import pandas as pd
import xarray as xr
import subprocess
import yaml
from ALLCools.utilities import genome_region_chunks
from ALLCools.mcds.utilities import write_ordered_chunks
from concurrent.futures import ProcessPoolExecutor, as_completed
from .rms_test import (
    permute_root_mean_square_test,
    calculate_residual,
    downsample_table,
)


def _read_single_allc(path, region):
    grange = region.split(":")[1]
    start, _ = grange.split("-")

    rows = []
    with pysam.TabixFile(path) as allc:
        for row in allc.fetch(region):
            _row = row.split("\t")[:6]
            rows.append(_row)

    data = pd.DataFrame(
        rows, columns=["chrom", "pos", "strand", "context", "mc", "cov"]
    ).set_index("pos")
    data.index = data.index.astype(np.int64)
    data["mc"] = data["mc"].astype(np.float64)
    data["cov"] = data["cov"].astype(np.float64)
    # important, turn the second column from cov to uc
    data["c"] = data["cov"] - data["mc"]

    # boarder case: if the pos equal to region start,
    # it might be included twice in adjacent regions
    # because both end and start are the same for single C
    data = data[data.index > int(start)]

    counts = data[["mc", "c"]].astype(np.float64)
    contexts = data["context"]
    return counts, contexts


def _make_count_table(path_list, region):
    all_contexts = {}
    all_counts = []
    for path in path_list:
        # counts has two columns: mc and c
        counts, contexts = _read_single_allc(path, region)
        all_contexts.update(contexts.to_dict())
        all_counts.append(counts)

    all_contexts = pd.Series(all_contexts, dtype="object").sort_index()
    total_counts = pd.concat(
        [df.reindex(all_contexts.index) for df in all_counts], axis=1
    ).fillna(0)
    return all_contexts, total_counts


def _perform_rms_batch(
    output_dir, allc_paths, samples, region, max_row_count, n_permute, min_pvalue
):
    contexts, total_counts = _make_count_table(path_list=allc_paths, region=region)

    n_samples = len(samples)
    total_values = {}
    p_values = {}

    if total_counts.shape[0] == 0:
        return None
    print(f"RMS tests for {total_counts.shape[0]} sites.")

    for pos, row in total_counts.iterrows():
        # table has two columns: mc and c
        table = row.values.reshape(len(samples), 2)
        _table = downsample_table(table, max_row_count=max_row_count)
        p = permute_root_mean_square_test(
            _table, n_permute=n_permute, min_pvalue=min_pvalue
        )
        # all the p-values will be saved for FDR
        p_values[pos] = p
        # but not all the values
        if p <= min_pvalue:
            residual = calculate_residual(table)
            # mc_residual = -c_residual, so only save one column
            residual = residual[:, [0]]
            table[:, 1] = table.sum(axis=1)  # turn c back to cov
            site_values = np.concatenate([table, residual], axis=1)
            total_values[pos] = site_values

    n_sites = len(total_values)
    site_matrix = np.zeros((n_samples, n_sites, 3))
    pos_order = []
    for i, (pos, site) in enumerate(total_values.items()):
        site_matrix[:, i, :] = site
        pos_order.append(pos)

    da = xr.DataArray(
        site_matrix[:, :, :2],
        coords=[samples, pos_order, ["mc", "cov"]],
        dims=["sample", "pos", "count_type"],
    )
    residual_da = xr.DataArray(
        site_matrix[:, :, 2], coords=[samples, pos_order], dims=["sample", "pos"]
    )
    if da.get_index("pos").size == 0:
        # pos dtype is wrong in case there is not record at all
        da.coords["pos"] = da.coords["pos"].astype(int)

    p_values = pd.Series(p_values, dtype="float64")
    p_values.index.name = "pos"
    da.coords["p-values"] = p_values

    contexts = pd.Series(contexts, dtype="object")
    contexts.index.name = "pos"
    da.coords["contexts"] = contexts

    da.coords["chrom"] = pd.Series(region.split(":")[0], index=da.get_index("pos"))
    ds = xr.Dataset({"dms_da": da, "dms_residual": residual_da})
    with np.errstate(divide="ignore"):
        ds["dms_da_frac"] = da.sel(count_type="mc") / da.sel(count_type="cov")

    # concatenate results
    # sort dataset and reset index
    ds = ds.sortby([ds.coords["chrom"], ds.coords["pos"]])
    ds = ds.rename({"pos": "dms"})
    # pos is still the genome coords
    ds.coords["pos"] = ds.coords["dms"].copy()
    # set str coords otherwise the zarr saving raise error
    ds.coords["chrom"] = ds.coords["chrom"].astype("str")
    # reset index to dms_id with int range
    ds.coords["dms"] = (
        ds.coords["chrom"].to_pandas().astype(str)
        + "-"
        + ds.coords["dms"].to_pandas().astype(str)
    )
    ds.coords["contexts"] = ds.coords["contexts"].astype("str")

    # rename none dimensional coords to prevent collision when merge with other ds
    ds = ds.rename({k: f"dms_{k}" for k in ds.coords.keys() if k not in ds.dims})

    # save total dms results
    output_path = f"{output_dir}/{region}.zarr"
    ds.to_zarr(output_path, mode="w")
    return output_path


def call_dms(
    output_dir,
    allc_paths,
    samples,
    chrom_size_path,
    cpu=1,
    max_row_count=50,
    n_permute=3000,
    min_pvalue=0.01,
    region=None,
):
    """
    Call DMS from multiple ALLC files

    Parameters
    ----------
    output_dir
    allc_paths
    samples
    chrom_size_path
    cpu
    max_row_count
    n_permute
    min_pvalue
    region

    Returns
    -------

    """
    pathlib.Path(output_dir).mkdir(exist_ok=True)

    with open(f"{output_dir}/.ALLCools", "w") as f:
        config = {"region_dim": "dms", "ds_region_dim": {"dms": "dms"}}
        yaml.dump(config, f)

    subprocess.run(
        f"cp {chrom_size_path} {output_dir}/chrom_sizes.txt", shell=True, check=True
    )
    chrom_size_path = f"{output_dir}/chrom_sizes.txt"

    # separate chrom chunks for parallel
    if region is None:
        regions = genome_region_chunks(chrom_size_path, bin_length=20000000)
    else:
        # only calculate user provided regions
        if isinstance(region, list):
            regions = region
        else:
            regions = [region]
            cpu = 1

    # temp dir
    dms_chunk_dir = pathlib.Path(f"{output_dir}/.dms_chunks")
    dms_chunk_dir.mkdir(exist_ok=True, parents=True)
    dms_dir = f"{output_dir}/dms"

    # trigger the numba JIT compilation before multiprocessing
    table = np.array([[0, 1], [0, 1]])
    permute_root_mean_square_test(table)
    calculate_residual(table)
    downsample_table(table, max_row_count=max_row_count)

    # parallel each chunk
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chunk_id, region in enumerate(regions):
            future = exe.submit(
                _perform_rms_batch,
                output_dir=dms_chunk_dir,
                allc_paths=allc_paths,
                samples=samples,
                region=region,
                max_row_count=max_row_count,
                n_permute=n_permute,
                min_pvalue=min_pvalue,
            )
            futures[future] = chunk_id

        chunks_to_write = {}
        for future in as_completed(futures):
            chunk_i = futures[future]
            output_path = future.result()
            chunks_to_write[chunk_i] = output_path

    write_ordered_chunks(
        chunks_to_write,
        final_path=dms_dir,
        append_dim="dms",
        engine="zarr",
        coord_dtypes=None,
    )

    subprocess.run(f"rm -rf {dms_chunk_dir}", shell=True)
    return
