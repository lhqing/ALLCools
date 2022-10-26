import pathlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import xarray as xr

from ..mcds import BaseDSChrom
from .rms_test import (
    calculate_residual,
    downsample_table,
    permute_root_mean_square_test,
    permute_root_mean_square_test_and_estimate_p,
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

    # sum the pos and neg
    # union the pos and neg index, otherwise if there is only one strand occur in the base_ds,
    # the other strand will be dropped by sum (xarray take the intersection of index)
    pos_index = pos_ds.get_index("pos")
    neg_index = neg_ds.get_index("pos")
    all_positions = pos_index.union(neg_index)
    cg_ds = pos_ds.reindex(pos=all_positions, fill_value=0) + neg_ds.reindex(pos=all_positions, fill_value=0)
    return cg_ds


def _core_dms_function(count_table, max_row_count, max_total_count, n_permute, min_pvalue, estimate_p=False):
    """Calculate dms p-value and residual for a single CpG count table."""
    count_table = downsample_table(count_table, max_row_count=max_row_count, max_total_count=max_total_count)
    if count_table.sum() == 0:
        return 1, np.zeros(count_table.shape[0])
    if estimate_p:
        p_value = permute_root_mean_square_test_and_estimate_p(
            table=count_table, n_permute=n_permute, min_pvalue=min_pvalue
        )
    else:
        p_value = permute_root_mean_square_test(table=count_table, n_permute=n_permute, min_pvalue=min_pvalue)

    residual = calculate_residual(count_table)
    # mc_residual = -c_residual, so only save one column
    residual = residual[:, 1]
    return p_value, residual


def add_groups_to_base_ds(base_ds, groups):
    """Add groups to base_ds.coords"""
    if isinstance(groups, (str, pathlib.Path)):
        groups = pd.read_csv(groups, header=None, index_col=0, names=["sample_id", "group"]).squeeze()

    # select samples
    ds_sample_id = base_ds.get_index("sample_id")
    use_sample = ds_sample_id.isin(groups.index)
    if use_sample.sum() != use_sample.size:
        # noinspection PyUnresolvedReferences
        base_ds = base_ds.sel(sample_id=use_sample)

    # add group info
    base_ds.coords["group"] = groups
    return base_ds


def call_dms_worker(
    groups,
    base_ds,
    mcg_pattern,
    n_permute,
    alpha,
    estimate_p,
    max_row_count,
    max_total_count,
    cpu,
    chrom,
    filter_sig,
    merge_strand,
    output_path,
    **output_kwargs,
):
    """
    Worker function for call_dms_from_base_ds.

    See call_dms_from_base_ds for parameter description.
    """
    if groups is not None:
        base_ds = add_groups_to_base_ds(base_ds, groups)

    # select CpG sites
    if mcg_pattern is not None:
        cg_base_ds = base_ds.select_mc_type(mcg_pattern)
        if merge_strand:
            # merge pos and neg strand
            cg_base_ds = _merge_pos_neg_cpg(cg_base_ds)
    else:
        cg_base_ds = base_ds

    if groups is not None:
        # group by group and sum base counts
        group_cg_base_ds = cg_base_ds[["data"]].groupby("group").sum(dim="sample_id").load()
    else:
        # if no group, load all data and call DMR across samples, the sample_id is renamed to group
        group_cg_base_ds = cg_base_ds[["data"]].rename({"sample_id": "group"}).load()

    pos_index = group_cg_base_ds.get_index("pos")
    data = group_cg_base_ds["data"]
    groups = group_cg_base_ds.get_index("group")
    with ProcessPoolExecutor(cpu) as exe:
        print(f"Perform DMS test for {pos_index.size} CpG sites and {groups.size} samples using {cpu} cores.")

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
                min_pvalue=alpha,
                estimate_p=estimate_p,
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
    residuals = pd.DataFrame.from_dict(residuals, orient="index", columns=groups)
    residuals.index.name = "pos"
    dms_ds = xr.Dataset({"dms_residual": residuals.astype("float32")})

    # add p-value
    p_values = pd.Series(p_values).astype("float32")
    p_values.index.name = "pos"
    p_values = p_values.reindex(pos_index)
    assert p_values.isna().sum() == 0, "Some positions are missing p-values."

    dms_ds.coords["p-values"] = p_values

    n_dms = dms_ds.get_index("pos").size

    if n_dms > 0:
        if estimate_p:
            # FDR correction, only valid for estimate_p
            from statsmodels.stats.multitest import multipletests

            _, q, *_ = multipletests(p_values, method="fdr_bh")
            q = pd.Series(q, index=p_values.index, name="q-values").astype("float32")
            dms_ds.coords["q-values"] = q

        # filter by p-value
        if filter_sig:
            if "q-values" in dms_ds.coords:
                dms_ds = dms_ds.sel(pos=dms_ds.coords["q-values"] <= alpha)
            else:
                dms_ds = dms_ds.sel(pos=dms_ds["p-values"] <= alpha)
    else:
        if estimate_p:
            empty_q = pd.Series([])
            empty_q.index.name = "pos"
            dms_ds.coords["q-values"] = empty_q

    if output_path is not None:
        dms_ds.to_zarr(output_path, **output_kwargs)
        return None
    else:
        return dms_ds


def call_dms_from_base_ds(
    base_ds_path,
    chrom,
    start,
    end,
    codebook_path=None,
    output_path=None,
    groups=None,
    mcg_pattern="CGN",
    cpu=1,
    n_permute=3000,
    alpha=0.01,
    max_row_count=50,
    max_total_count=3000,
    filter_sig=True,
    merge_strand=True,
    estimate_p=True,
    **output_kwargs,
):
    """
    Call DMS for a single genome region.

    Parameters
    ----------
    base_ds_path :
        Path to the BaseDS, wildcard accepted if data stored in multiple BaseDS.
    chrom :
        Chromosome name.
    start :
        Start position.
    end :
        End position.
    codebook_path :
        Path to the mc_type codebook dataset.
    output_path :
        Path to the output DMS dataset.
        If provided, the result will be saved to disk.
        If not, the result will be returned.
    groups :
        Grouping information for the samples.
        If None, perform DMS test on all samples in the BaseDS.
        If provided, first group the samples by the group information, then perform DMS test on each group.
        Samples not occur in the group information will be ignored.
    mcg_pattern :
        Pattern of the methylated cytosine, default is "CGN".
    cpu :
        Number of CPU to use.
    n_permute :
        Number of permutation to perform.
    alpha :
        Minimum p-value/q-value to consider a site as significant.
    max_row_count :
        Maximum number of base counts for each row (sample) in the DMS input count table.
    max_total_count :
        Maximum total number of base counts in the DMS input count table.
    filter_sig :
        Whether to filter out the non-significant sites in output DMS dataset.
    merge_strand :
            Whether to merge the base counts of CpG sites next to each other.
    estimate_p :
        Whether to estimate p-value by approximate the null distribution of S as normal distribution.
        The resolution of the estimated p-value is much higher than the exact p-value,
        which is necessary for multiple testing correction.
        FDR corrected q-value is also estimated if this option is enabled.
    output_kwargs :
        Keyword arguments for the output DMS dataset, pass to xarray.Dataset.to_zarr.

    Returns
    -------
    xarray.Dataset if output_path is None, otherwise None.
    """
    base_ds = BaseDSChrom.open(
        f"{base_ds_path}/{chrom}",
        codebook_path=f"{codebook_path}/{chrom}",
        start=start,
        end=end,
    )

    dms_ds = call_dms_worker(
        groups=groups,
        base_ds=base_ds,
        mcg_pattern=mcg_pattern,
        n_permute=n_permute,
        alpha=alpha,
        estimate_p=estimate_p,
        max_row_count=max_row_count,
        max_total_count=max_total_count,
        cpu=cpu,
        chrom=chrom,
        filter_sig=filter_sig,
        merge_strand=merge_strand,
        output_path=output_path,
        **output_kwargs,
    )
    return dms_ds
