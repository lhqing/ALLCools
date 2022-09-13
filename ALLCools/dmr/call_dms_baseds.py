from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import xarray as xr

from ..mcds import BaseDSChrom
from .rms_test import (
    calculate_residual,
    downsample_table,
    permute_root_mean_square_test,
)


def _merge_pos_neg_cpg(cg_ds):
    # assume the ds already gone through CpG selection, only CG sites remained
    strand = cg_ds["codebook"].sum(dim="mc_type").to_pandas()
    neg_pos = strand < 0

    # separate pos and neg
    # noinspection PyUnresolvedReferences
    pos_ds = cg_ds.sel(pos=(~neg_pos).values)
    neg_ds = cg_ds.sel(pos=neg_pos.values)

    # modify the pos by -1 so the neg_ds and pos_ds has the same pos,
    # and can be sum together
    neg_ds["pos"] = neg_ds["pos"].to_pandas() - 1
    cg_ds = pos_ds + neg_ds
    return cg_ds


def _core_dms_function(count_table, max_row_count, max_total_count, n_permute, min_pvalue):
    """Calculate dms p-value and residual for a single CpG count table."""
    count_table = downsample_table(count_table, max_row_count=max_row_count, max_total_count=max_total_count)
    p_value = permute_root_mean_square_test(table=count_table, n_permute=n_permute, min_pvalue=min_pvalue)
    residual = calculate_residual(count_table)
    # mc_residual = -c_residual, so only save one column
    residual = residual[:, 1]
    return p_value, residual


def _call_dms_worker(
    base_ds_path,
    codebook_path,
    output_path,
    chrom,
    start,
    end,
    groups_path=None,
    mcg_pattern="CGN",
    cpu=1,
    n_permute=3000,
    min_pvalue=0.034,
    max_row_count=50,
    max_total_count=3000,
    **output_kwargs,
):
    base_ds = BaseDSChrom.open(
        f"{base_ds_path}/{chrom}",
        codebook_path=f"{codebook_path}/{chrom}",
        start=start,
        end=end,
    )

    if groups_path is not None:
        sample_group = pd.read_csv(groups_path, header=None, index_col=0, names=["sample_id", "group"]).squeeze()
        base_ds.coords["group"] = sample_group

    # select CpG sites
    cg_base_ds = base_ds.select_mc_type(mcg_pattern)

    # merge pos and neg strand
    cg_base_ds = _merge_pos_neg_cpg(cg_base_ds)

    if groups_path is not None:
        # group by group and sum base counts
        group_cg_base_ds = cg_base_ds[["data"]].groupby("group").sum(dim="sample_id").load()
    else:
        # if no group, load all data and call DMR across samples, the sample_id is renamed to group
        group_cg_base_ds = cg_base_ds[["data"]].rename({"sample_id": "group"}).load()

    with ProcessPoolExecutor(cpu) as exe:
        pos_index = group_cg_base_ds.get_index("pos")
        data = group_cg_base_ds["data"]

        futures = {}
        for pos in pos_index:
            count_table = data.sel(pos=pos).to_pandas()
            count_table["uc"] = count_table["cov"] - count_table["mc"]
            count_table = count_table[["mc", "uc"]].values

            future = exe.submit(
                _core_dms_function,
                count_table=count_table,
                max_row_count=max_row_count,
                max_total_count=max_total_count,
                n_permute=n_permute,
                min_pvalue=min_pvalue,
            )
            futures[future] = pos

        p_values = {}
        residuals = {}
        for future in as_completed(futures):
            pos = futures[future]
            try:
                p_value, residual = future.result()
            except Exception as e:
                print(f"Position {chrom} {pos} got an error.")
                raise e

            p_values[pos] = p_value
            residuals[pos] = residual

    # create ds
    residuals = pd.DataFrame.from_dict(residuals, orient="index", columns=group_cg_base_ds.get_index("group"))
    residuals.index.name = "pos"
    dms_ds = xr.Dataset({"dms_residual": residuals})

    # add p-value
    p_values = pd.Series(p_values)
    dms_ds.coords["p-values"] = pd.Series(p_values, index=pos_index)

    # filter by p-value
    dms = dms_ds.sel(pos=dms_ds["p-values"] < 0.1)
    dms.to_zarr(output_path, **output_kwargs)
    return
