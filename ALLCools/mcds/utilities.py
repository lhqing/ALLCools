import glob
import json
import logging
import pathlib
import subprocess
import warnings
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr
import yaml

log = logging.getLogger()


def calculate_posterior_mc_frac(mc_da, cov_da, var_dim=None, normalize_per_cell=True, clip_norm_value=10):
    # so we can do post_frac only in a very small set of gene to prevent memory issue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # here we expected to see a true_divide warning due to cov=0
        raw_frac = mc_da / cov_da

    if isinstance(raw_frac, np.ndarray):
        # np.ndarray
        ndarray = True
    else:
        ndarray = False

    if ndarray:
        cell_rate_mean = np.nanmean(raw_frac, axis=1)
        cell_rate_var = np.nanvar(raw_frac, axis=1)
    else:
        # assume xr.DataArray
        if var_dim is None:
            cell_rate_mean = raw_frac.mean(axis=1)  # this skip na
            cell_rate_var = raw_frac.var(axis=1)  # this skip na
        else:
            cell_rate_mean = raw_frac.mean(dim=var_dim)  # this skip na
            cell_rate_var = raw_frac.var(dim=var_dim)  # this skip na

    # based on beta distribution mean, var
    # a / (a + b) = cell_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = cell_rate_var
    # calculate alpha beta value for each cell
    cell_a = (1 - cell_rate_mean) * (cell_rate_mean**2) / cell_rate_var - cell_rate_mean
    cell_b = cell_a * (1 / cell_rate_mean - 1)

    # cell specific posterior rate
    post_frac: Union[np.ndarray, xr.DataArray]
    if ndarray:
        post_frac = (mc_da + cell_a[:, None]) / (cov_da + cell_a[:, None] + cell_b[:, None])
    else:
        post_frac = (mc_da + cell_a) / (cov_da + cell_a + cell_b)

    if normalize_per_cell:
        # there are two ways of normalizing per cell, by posterior or prior mean:
        # prior_mean = cell_a / (cell_a + cell_b)
        # posterior_mean = post_rate.mean(dim=var_dim)

        # Here I choose to use prior_mean to normalize cell,
        # therefore all cov == 0 features will have normalized rate == 1 in all cells.
        # i.e. 0 cov feature will provide no info
        prior_mean = cell_a / (cell_a + cell_b)
        if ndarray:
            post_frac = post_frac / prior_mean[:, None]
        else:
            post_frac = post_frac / prior_mean
        if clip_norm_value is not None:
            if isinstance(post_frac, np.ndarray):
                # np.ndarray
                post_frac[post_frac > clip_norm_value] = clip_norm_value
            else:
                # xarray.DataArray
                post_frac = post_frac.where(post_frac < clip_norm_value, clip_norm_value)
    return post_frac


def calculate_posterior_mc_frac_lazy(
    mc_da,
    cov_da,
    var_dim,
    output_prefix,
    cell_chunk=20000,
    dask_cell_chunk=500,
    normalize_per_cell=True,
    clip_norm_value=10,
):
    """
    Calculate posterior mc frac for each cell, dask version.

    Running calculate_posterior_mc_rate with dask array and directly save to disk.
    This is highly memory efficient. Use this for dataset larger than machine memory.
    """
    cell_list = mc_da.get_index("cell")
    cell_chunks = [cell_list[i : i + cell_chunk] for i in range(0, cell_list.size, cell_chunk)]

    output_paths = []
    for chunk_id, cell_list_chunk in enumerate(cell_chunks):
        _mc_da = mc_da.sel(cell=cell_list_chunk)
        _cov_da = cov_da.sel(cell=cell_list_chunk)
        post_rate = calculate_posterior_mc_frac(
            mc_da=_mc_da,
            cov_da=_cov_da,
            var_dim=var_dim,
            normalize_per_cell=normalize_per_cell,
            clip_norm_value=clip_norm_value,
        )
        if len(cell_chunks) == 1:
            chunk_id = ""
        else:
            chunk_id = f".{chunk_id}"
        output_path = output_prefix + f".{var_dim}_da_frac{chunk_id}.mcds"

        # to_netcdf trigger the dask computation, and save output directly into disk, quite memory efficient
        post_rate.to_netcdf(output_path)
        output_paths.append(output_path)

    chunks = {"cell": dask_cell_chunk}
    total_post_rate = xr.concat([xr.open_dataarray(path, chunks=chunks) for path in output_paths], dim="cell")
    return total_post_rate


def calculate_gch_rate(mcds, var_dim="chrom100k"):
    rate_da = mcds.sel(mc_type=["GCHN", "HCHN"]).add_mc_rate(
        dim=var_dim, da=f"{var_dim}_da", normalize_per_cell=False, inplace=False
    )
    # (PCG - PCH) / (1 - PCH)
    real_gc_rate = (rate_da.sel(mc_type="GCHN") - rate_da.sel(mc_type="HCHN")) / (1 - rate_da.sel(mc_type="HCHN"))
    real_gc_rate = real_gc_rate.transpose("cell", var_dim).values
    real_gc_rate[real_gc_rate < 0] = 0

    # norm per cell
    cell_overall_count = mcds[f"{var_dim}_da"].sel(mc_type=["GCHN", "HCHN"]).sum(dim=var_dim)
    cell_overall_rate = cell_overall_count.sel(count_type="mc") / cell_overall_count.sel(count_type="cov")
    gchn = cell_overall_rate.sel(mc_type="GCHN")
    hchn = cell_overall_rate.sel(mc_type="HCHN")
    overall_gchn = (gchn - hchn) / (1 - hchn)
    real_gc_rate = real_gc_rate / overall_gchn.values[:, None]
    return real_gc_rate


def get_mean_dispersion(x, obs_dim):
    # mean
    mean = x.mean(dim=obs_dim).load()
    abs_mean = np.abs(mean)

    # var
    # enforce R convention (unbiased estimator) for variance
    mean_sq = (x * x).mean(dim=obs_dim).load()
    var = (mean_sq - abs_mean**2) * (x.sizes[obs_dim] / (x.sizes[obs_dim] - 1))

    # now actually compute the dispersion
    abs_mean.where(abs_mean > 1e-12, 1e-12)  # set entries equal to zero to small value
    # raw dispersion is the variance normalized by mean
    dispersion = var / abs_mean
    return mean, dispersion


def highly_variable_methylation_feature(
    cell_by_feature_matrix,
    feature_mean_cov,
    obs_dim=None,
    var_dim=None,
    min_disp=0.5,
    max_disp=None,
    min_mean=0,
    max_mean=5,
    n_top_feature=None,
    bin_min_features=5,
    mean_binsize=0.05,
    cov_binsize=100,
):
    """
    Calculate highly variable methylation features.

    Adapted from Scanpy, the main difference is that,
    this function normalize dispersion based on both mean and cov bins.
    """
    # RNA is only scaled by mean, but methylation is scaled by both mean and cov
    log.info("extracting highly variable features")

    if n_top_feature is not None:
        log.info(
            "If you pass `n_top_feature`, all cutoffs are ignored. " "Features are ordered by normalized dispersion."
        )

    # warning for extremely low cov
    low_cov_portion = (feature_mean_cov < 10).sum() / feature_mean_cov.size
    if low_cov_portion > 0.2:
        log.warning(
            f"{int(low_cov_portion * 100)}% feature with < 10 mean cov, "
            f"consider filter by cov before find highly variable feature. "
            f"Otherwise some low coverage feature may be elevated after normalization."
        )

    if len(cell_by_feature_matrix.dims) != 2:
        raise ValueError(
            f"Input cell_by_feature_matrix is not 2-D matrix, " f"got {len(cell_by_feature_matrix.dims)} dim(s)"
        )
    else:
        if (obs_dim is None) or (var_dim is None):
            obs_dim, var_dim = cell_by_feature_matrix.dims

    # rename variable
    x = cell_by_feature_matrix
    cov = feature_mean_cov

    mean, dispersion = get_mean_dispersion(x, obs_dim=obs_dim)
    dispersion = np.log(dispersion + 1e-12)

    # all the following quantities are "per-feature" here
    df = pd.DataFrame(index=cell_by_feature_matrix.get_index(var_dim))
    df["mean"] = mean.to_pandas().copy()
    df["dispersion"] = dispersion.to_pandas().copy()
    df["cov"] = cov.to_pandas().copy()

    # instead of n_bins, use bin_size, because cov and mc are in different scale
    df["mean_bin"] = (df["mean"] / mean_binsize).astype(int)
    df["cov_bin"] = (df["cov"] / cov_binsize).astype(int)

    # save bin_count df, gather bins with more than bin_min_features features
    bin_count = (
        df.groupby(["mean_bin", "cov_bin"]).apply(lambda i: i.shape[0]).reset_index().sort_values(0, ascending=False)
    )
    bin_count.head()
    bin_more_than = bin_count[bin_count[0] > bin_min_features]
    if bin_more_than.shape[0] == 0:
        raise ValueError(f"No bin have more than {bin_min_features} features, use larger bin size.")

    # for those bin have too less features, merge them with closest bin in manhattan distance
    # this usually don't cause much difference (a few hundred features), but the scatter plot will look more nature
    index_map = {}
    for _, (mean_id, cov_id, count) in bin_count.iterrows():
        if count > 1:
            index_map[(mean_id, cov_id)] = (mean_id, cov_id)
        manhattan_dist = (bin_more_than["mean_bin"] - mean_id).abs() + (bin_more_than["cov_bin"] - cov_id).abs()
        closest_more_than = manhattan_dist.sort_values().index[0]
        closest = bin_more_than.loc[closest_more_than]
        index_map[(mean_id, cov_id)] = tuple(closest.tolist()[:2])
    # apply index_map to original df
    raw_bin = df[["mean_bin", "cov_bin"]].set_index(["mean_bin", "cov_bin"])
    raw_bin["use_mean"] = pd.Series(index_map).apply(lambda i: i[0])
    raw_bin["use_cov"] = pd.Series(index_map).apply(lambda i: i[1])
    df["mean_bin"] = raw_bin["use_mean"].values
    df["cov_bin"] = raw_bin["use_cov"].values

    # calculate bin mean and std, now disp_std_bin shouldn't have NAs
    disp_grouped = df.groupby(["mean_bin", "cov_bin"])["dispersion"]
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)

    # actually do the normalization
    _mean_norm = disp_mean_bin.loc[list(zip(df["mean_bin"], df["cov_bin"]))]
    _std_norm = disp_std_bin.loc[list(zip(df["mean_bin"], df["cov_bin"]))]
    df["dispersion_norm"] = (
        df["dispersion"].values - _mean_norm.values  # use values here as index differs
    ) / _std_norm.values
    dispersion_norm = df["dispersion_norm"].values.astype("float32")

    # Select n_top_feature
    if n_top_feature is not None:
        feature_subset = df.index.isin(df.sort_values("dispersion_norm", ascending=False).index[:5000])
    else:
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        feature_subset = np.logical_and.reduce(
            (
                mean > min_mean,
                mean < max_mean,
                dispersion_norm > min_disp,
                dispersion_norm < max_disp,
            )
        )
    df["feature_select"] = feature_subset
    log.info("    finished")
    return df


def determine_engine(dataset_paths):
    def _single_path(path):
        if pathlib.Path(f"{path}/.zgroup").exists():
            e = "zarr"
        else:
            e = None  # default for None is netcdf4 or xarray will guess
        return e

    def _multi_paths(paths):
        engines = []
        for path in paths:
            e = _single_path(path)
            engines.append(e)
        _engine = list(set(engines))
        if len(_engine) > 1:
            raise ValueError(
                f'Can not open a mixture of netcdf4 and zarr files in "{dataset_paths}", please use '
                f"`allcools convert-mcds-to-zarr {{dataset_paths}}` to convert the storage of all MCDS files "
                f"to zarr or netcdf. We recommend you use the zarr format, which gives better IO "
                f"performance according to our experience."
            )
        else:
            _engine = _engine[0]
        return _engine

    if isinstance(dataset_paths, (str, pathlib.PosixPath)):
        if "*" in str(dataset_paths):
            engine = _multi_paths(list(glob.glob(dataset_paths)))
        else:
            # single mcds path
            engine = _single_path(dataset_paths)
    else:
        if len(dataset_paths) == 1:
            engine = _single_path(dataset_paths[0])
        else:
            engine = _multi_paths(dataset_paths)
    return engine


def obj_to_str(ds, coord_dtypes=None):
    if coord_dtypes is None:
        coord_dtypes = {}

    for k, v in ds.coords.items():
        if np.issubdtype(v, np.object) or np.issubdtype(v, np.unicode):
            if k in coord_dtypes:
                ds.coords[k] = v.load().astype(coord_dtypes[k])
            else:
                ds.coords[k] = v.load().astype(str)
    return


def write_ordered_chunks(
    chunks_to_write,
    final_path,
    append_dim,
    engine="zarr",
    coord_dtypes=None,
    dtype=None,
):
    # some function may return None if the chunk is empty
    chunks_to_write = {k: v for k, v in chunks_to_write.items() if v is not None}

    # open all chunks to determine string coords dtype lengths
    total_ds = xr.open_mfdataset(
        list(chunks_to_write.values()),
        concat_dim=append_dim,
        combine="nested",
        engine=engine,
    )
    obj_to_str(total_ds, coord_dtypes)
    # string dtype is set to coord maximum length, so need to load all coords to determine
    # otherwise the chunk writing will truncate string coords if the dtype length is wrong
    coord_dtypes = {k: v.dtype for k, v in total_ds.coords.items()}
    del total_ds

    # write chunks
    wrote = False
    for _, output_path in sorted(chunks_to_write.items(), key=lambda _i: _i[0]):
        # save chunk into the zarr
        chunk_ds = xr.open_dataset(output_path, engine=engine).load()
        obj_to_str(chunk_ds, coord_dtypes)
        if dtype is not None:
            chunk_ds = chunk_ds.astype(dtype)
        if not wrote:
            wrote = True
            # create the new da
            chunk_ds.to_zarr(final_path, mode="w")
        else:
            chunk_ds.to_zarr(final_path, mode="a", append_dim=append_dim)
    return


def convert_to_zarr(paths):
    """Convert xarray.Dataset stored in other backends into zarr backend."""

    def _convert_single_path(p):
        if not pathlib.Path(p).exists():
            raise FileNotFoundError(f"{p} not exist.")

        tmp_p = f"{p}_convert_tmp"
        if determine_engine(tmp_p) != "zarr":
            ds = xr.open_dataset(p)
            # this will load the whole dataset
            ds.to_zarr(tmp_p)
            try:
                subprocess.run(["mv", p, f"{p}_to_delete"], check=True, stderr=subprocess.PIPE, encoding="utf8")
                subprocess.run(["mv", tmp_p, p], check=True, stderr=subprocess.PIPE, encoding="utf8")
                subprocess.run(["rm", "-rf", f"{p}_to_delete"], check=True, stderr=subprocess.PIPE, encoding="utf8")
            except subprocess.CalledProcessError as e:
                print(e.stderr)
                raise e
        return

    if isinstance(paths, (str, pathlib.PosixPath)):
        paths = str(paths)
        if "*" in paths:
            for path in glob.glob(paths):
                _convert_single_path(path)
        else:
            _convert_single_path(paths)
    else:
        for path in paths:
            convert_to_zarr(path)
    return


def update_dataset_config(
    output_dir, add_ds_region_dim=None, change_region_dim=None, config=None, add_ds_sample_dim=None
):
    # update RegionDS default dimension
    try:
        with open(f"{output_dir}/.ALLCools") as f:
            _config = yaml.load(f, yaml.SafeLoader)
    except FileNotFoundError:
        _config = {"region_dim": None, "ds_region_dim": {}, "ds_sample_dim": {}}

    if config is not None:
        _config.update(config)

    if change_region_dim is not None:
        _config["region_dim"] = change_region_dim

    if add_ds_region_dim is not None:
        _config["ds_region_dim"].update(add_ds_region_dim)

    if add_ds_sample_dim is not None:
        _config["ds_sample_dim"].update(add_ds_sample_dim)

    with open(f"{output_dir}/.ALLCools", "w") as f:
        yaml.dump(_config, f)
    return


def reduce_zarr_coords_chunks(ds_path, max_size=10000000):
    """Reduce zarr dataset coordinates chunks to speed up coordinate loading."""
    ds = xr.open_zarr(ds_path)

    # update coords encoding
    new_coords_encoding = {}
    for k, v in ds.coords.items():
        encoding = v.encoding.copy()
        chunks = []
        preferred_chunks = {}
        for dim in v.dims:
            chunk_size = min(ds.get_index(dim).size, max_size)
            chunks.append(chunk_size)
            preferred_chunks[dim] = chunk_size
        encoding["chunks"] = tuple(chunks)
        encoding["preferred_chunks"] = preferred_chunks
        try:
            encoding["compressor"].clevel = 1
        except KeyError:
            pass
        new_coords_encoding[k] = encoding

    for k, encoding in new_coords_encoding.items():
        ds.coords[k].encoding = encoding

    # create a coords only dataset, save to temp zarr
    coord_ds = xr.Dataset({}, coords=ds.coords).load()
    # turn obj to str
    for k, v in coord_ds.coords.items():
        if v.dtype == "O":
            coord_ds.coords[k] = v.astype(str)
    temp_zarr_path = f"{ds_path}_temp.zarr"
    coord_ds.to_zarr(temp_zarr_path, mode="w")

    # move previous coords and zmetadata
    for coord in coord_ds.coords.keys():
        subprocess.run(["cp", "-r", f"{ds_path}/{coord}", f"{temp_zarr_path}/{coord}_backup"], check=True)
        subprocess.run(["rm", "-rf", f"{ds_path}/{coord}"], check=True)
    # move new coords
    for coord in coord_ds.coords.keys():
        subprocess.run(["cp", "-r", f"{temp_zarr_path}/{coord}", f"{ds_path}/{coord}"], check=True)

    # zmetadata path
    old_zmetadata_path = f"{ds_path}/.zmetadata"
    new_zmetadata_path = f"{temp_zarr_path}/.zmetadata"
    # load both zmetadata
    with open(old_zmetadata_path) as f:
        zmeta_tmp = json.load(f)
    with open(new_zmetadata_path) as f:
        zmeta_coords = json.load(f)
    # update coord zmetadata
    z_format = zmeta_tmp["zarr_consolidated_format"]
    if z_format == 1:
        update_meta = {k: v for k, v in zmeta_coords["metadata"].items() if not k.startswith(".")}
        for k in zmeta_tmp["metadata"].keys():
            if k in update_meta:
                zmeta_tmp["metadata"][k] = update_meta[k]
    else:
        raise NotImplementedError(f"parser of zarr_consolidated_format {z_format} not implemented")

    # save the changes
    with open(old_zmetadata_path, "w") as f:
        json.dump(zmeta_tmp, f)

    # delete temp zarr
    subprocess.run(["rm", "-rf", temp_zarr_path], check=True)
    return


def binnify(chromsizes, binsize):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    bins : :py:class:`pandas.DataFrame`
        Dataframe with columns: ``chrom``, ``start``, ``end``.
    """

    def _each(chrom):
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins + 1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame(
            {"chrom": [chrom] * n_bins, "start": binedges[:-1], "end": binedges[1:]},
            columns=["chrom", "start", "end"],
        )

    bintable = pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)

    bintable["chrom"] = pd.Categorical(bintable["chrom"], categories=list(chromsizes.index), ordered=True)

    return bintable
