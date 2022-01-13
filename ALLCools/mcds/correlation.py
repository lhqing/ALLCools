"""
See here https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes

Calculate correlation between two matrix, row by row
"""

from numba import njit
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from concurrent.futures import ProcessPoolExecutor, as_completed
import anndata
import scanpy as sc
import xarray as xr
from sklearn.impute import SimpleImputer


@njit
def _mean(a):
    n = len(a)
    b = np.empty(n)
    for i in range(n):
        b[i] = a[i].mean()
    return b


@njit
def _std(a):
    n = len(a)
    b = np.empty(n)
    for i in range(n):
        b[i] = a[i].std()
    return b


@njit
def _corr(a, b, row, col):
    """
    Correlation between rows in a and b, no nan value
    """
    _, k = a.shape

    mu_a = _mean(a)
    mu_b = _mean(b)
    sig_a = _std(a)
    sig_b = _std(b)

    out = np.zeros(shape=row.shape, dtype=np.float32)

    for idx in range(out.size):
        i = row[idx]
        j = col[idx]

        _sig_a = sig_a[i]
        _sig_b = sig_b[j]
        if _sig_a == 0 or _sig_b == 0:
            # if any variable std == 0
            out[idx] = np.nan
        else:
            out[idx] = (a[i] - mu_a[i]) @ (b[j] - mu_b[j]) / k / _sig_a / _sig_b
    return out


def _corr_preprocess(da, sample_mch, sample_mcg, cpu=1):
    imputer = SimpleImputer(copy=False)
    df = da.to_pandas()
    imputer.fit_transform(df)
    assert df.isna().values.sum() == 0

    adata = anndata.AnnData(
        X=df.values,
        obs=pd.DataFrame([], index=df.index),
        var=pd.DataFrame([], index=df.columns),
    )

    adata.var["pos"] = da["_pos"]

    # standardize global
    sample_mcg = (sample_mcg - sample_mcg.mean()) / sample_mcg.std()
    sample_mch = (sample_mch - sample_mch.mean()) / sample_mch.std()
    # add global
    adata.obs["sample_mch"] = sample_mch
    adata.obs["sample_mcg"] = sample_mcg

    # remove dmr that has 0 std
    # otherwise regress out will fail.
    adata = adata[:, adata.X.std(axis=0) > 0].copy()

    if adata.shape[1] > 0:
        # regress out global mCH and mCG
        # happens on each regions
        sc.pp.regress_out(adata, keys=["sample_mch", "sample_mcg"], n_jobs=cpu)
        sc.pp.scale(adata)
    else:
        print('All features are removed due to 0 std')
    return adata


def corr(
    adata_a,
    adata_b,
    max_dist,
    cpu=1,
    method="pearson",
    shuffle_sample=None,
    calculate_n=None,
):
    _adata_a = adata_a.copy()
    _adata_b = adata_b.copy()

    # trigger njit
    _corr(np.array([[0, 1]]), np.array([[0, 1]]), np.array([0]), np.array([0]))

    # this is used when null type is sample
    if shuffle_sample:
        obs_names = _adata_a.obs_names.values.copy()
        np.random.shuffle(obs_names)
        _adata_a.obs_names = obs_names
        _adata_a = _adata_a[_adata_b.obs_names, :]

    # create mask
    mask = coo_matrix(
        np.abs(
            (_adata_a.var["pos"].values[:, None] - _adata_b.var["pos"].values[None, :])
        )
        < max_dist
    )

    # this is used when calculating null, downsample corr call to save time
    if (calculate_n is not None) and (calculate_n < mask.data.size):
        rand_loc = np.random.choice(
            range(mask.data.size), size=int(calculate_n), replace=False
        )
        rand_loc = sorted(rand_loc)
        mask = coo_matrix(
            (mask.data[rand_loc], (mask.row[rand_loc], mask.col[rand_loc])),
            shape=mask.shape,
        )

    # adata is sample by region
    # _corr input is region by sample
    a = _adata_a.X.T
    b = _adata_b.X.T

    if method.lower()[0] == "p":
        pass
    elif method.lower()[0] == "s":
        # turn a, b in to rank matrix
        a = a.argsort(axis=1).argsort(axis=1)
        b = b.argsort(axis=1).argsort(axis=1)
    else:
        raise ValueError("Method can only be pearson or spearman")

    a = a.astype(np.float32)
    b = b.astype(np.float32)
    chunk_size = max(100, mask.data.size // cpu)

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for i, chunk_start in enumerate(range(0, mask.data.size, chunk_size)):
            future = exe.submit(
                _corr,
                a=a,
                b=b,
                row=mask.row[chunk_start : chunk_start + chunk_size],
                col=mask.col[chunk_start : chunk_start + chunk_size],
            )
            futures[future] = i

        total_data = {}
        for future in as_completed(futures):
            chunk_id = futures[future]
            total_data[chunk_id] = future.result()

    if len(total_data) != 0:
        mask.data = np.concatenate([total_data[k] for k in sorted(total_data.keys())])
    result = anndata.AnnData(
        X=mask.astype(np.float32).tocsr(), obs=_adata_a.var, var=_adata_b.var
    )
    return result


def _verify_correlation_da(data: xr.DataArray, region_dim, pos):
    if len(data.dims) != 2:
        raise ValueError(f"Only support 2D data array, got {len(data.dims)} dims")
    if region_dim is None:
        region_dim = data.dims[1]
        print(f"Using {region_dim} as region_dim")
    data.transpose(*(v for v in data.dims if v != region_dim), region_dim)

    # test if {region_dim}_chrom exist
    try:
        data.coords[f"{region_dim}_chrom"]
    except KeyError as e:
        raise e

    if pos is None:
        print(f"Using {region_dim}_start and {region_dim}_end to calculate pos")
        start = data.coords[f"{region_dim}_start"]
        end = data.coords[f"{region_dim}_end"]
        center = (end + start) // 2
        data.coords["_pos"] = center
    elif isinstance(pos, str):
        data.coords["_pos"] = data.coords[pos]
    else:
        data.coords["_pos"] = pos

    return data, region_dim


def region_correlation(
    data_a,
    data_b,
    sample_mch,
    sample_mcg,
    region_dim_a=None,
    region_dim_b=None,
    pos_a=None,
    pos_b=None,
    method="pearson",
    max_dist=1000000,
    cpu=1,
    null="sample",
    null_n=100000,
    chroms=None,
):
    data_a, region_dim_a = _verify_correlation_da(data_a, region_dim_a, pos_a)
    data_b, region_dim_b = _verify_correlation_da(data_b, region_dim_b, pos_b)

    total_results = []
    null_results = []
    region_null_pairs = 20

    for chrom_a, chrom_da_a in data_a.groupby(f"{region_dim_a}_chrom"):
        if (chroms is not None) and (chrom_a not in chroms):
            continue
        adata_a = _corr_preprocess(
            da=chrom_da_a, sample_mch=sample_mch, sample_mcg=sample_mcg, cpu=cpu
        )
        for chrom_b, chrom_da_b in data_b.groupby(f"{region_dim_b}_chrom"):
            if (chroms is not None) and (chrom_b not in chroms):
                continue
            adata_b = _corr_preprocess(
                da=chrom_da_b, sample_mch=sample_mch, sample_mcg=sample_mcg, cpu=cpu
            )

            # we use trans chroms to calculate region null
            if chrom_a != chrom_b:
                if (null == "region") & (region_null_pairs > 0):
                    print(f"Calculating region null {chrom_a} {chrom_b}")
                    null_corr_matrix = corr(
                        adata_a,
                        adata_b,
                        max_dist=999999999999,
                        method=method,
                        shuffle_sample=False,
                        calculate_n=null_n,
                        cpu=cpu,
                    )
                    null_results.append(null_corr_matrix.X.data)
                    region_null_pairs -= 1
                continue
            else:
                print(f"Calculating {chrom_a}")
                # this is calculating real corr
                corr_matrix = corr(
                    adata_a,
                    adata_b,
                    max_dist=max_dist,
                    method=method,
                    shuffle_sample=False,
                    calculate_n=None,
                    cpu=cpu,
                )
                total_results.append(corr_matrix)

                if null == "sample":
                    # this is calculating null corr
                    null_corr_matrix = corr(
                        adata_a,
                        adata_b,
                        max_dist=max_dist,
                        method=method,
                        shuffle_sample=True,
                        calculate_n=null_n,
                        cpu=cpu,
                    )
                    # no need to save index for null, just use its value
                    null_results.append(null_corr_matrix.X.data)
    true_results = anndata.concat(total_results, join="outer")
    # anndata somehow do not concat var and the var order is also changed
    true_results.var = pd.concat([adata.var for adata in total_results]).loc[
        true_results.var_names
    ]
    null_results = np.concatenate(null_results)
    return true_results, null_results


def _select_corr_cutoff(true: np.ndarray, null: np.ndarray, alpha=0.05, direction="+"):
    if direction == "both":
        return (
            _select_corr_cutoff(true, null, alpha / 2, direction="+"),
            _select_corr_cutoff(true, null, alpha / 2, direction="-"),
        )
    if direction == "+":
        cur_corr = 0.8
        corr_range = np.arange(0.8, 0.2, -0.01)
    else:
        cur_corr = -0.8
        corr_range = np.arange(-0.8, -0.2, 0.01)
    for corr_cutoff in corr_range:
        if direction == "+":
            # positive corr
            # p value of this corr based on null distribution
            p = (null > corr_cutoff).sum() / null.size
            # ture data having smaller P than this point
            p_smaller = (true > corr_cutoff).sum()
        else:
            # negative corr
            p = (null < corr_cutoff).sum() / null.size
            p_smaller = (true < corr_cutoff).sum()

        if p_smaller == 0:
            q = 0
        else:
            # q value with BH correction
            q = p / p_smaller * true.size
        # print(corr_cutoff, q)
        if q > alpha:
            break
        cur_corr = corr_cutoff
    return cur_corr


def get_corr_table(
    total_results, null_results, region_dim_a, region_dim_b, direction="-", alpha=0.05
):
    true = total_results.X.data
    null = null_results
    total_results.X = total_results.X.tocoo()

    if direction == "+":
        cutoff = _select_corr_cutoff(true, null, alpha=alpha, direction="+")
        selected = total_results.X.data > cutoff
        print(f"Correlation cutoff corr > {cutoff}")
    elif direction == "-":
        cutoff = _select_corr_cutoff(true, null, alpha=alpha, direction="-")
        selected = total_results.X.data < cutoff
        print(f"Correlation cutoff corr < {cutoff}")
    elif direction == "both":
        pos_cutoff, neg_cutoff = _select_corr_cutoff(
            true, null, alpha=alpha, direction="both"
        )
        selected = (total_results.X.data > pos_cutoff) | (
            total_results.X.data < neg_cutoff
        )
        print(f"Correlation cutoff (corr > {pos_cutoff}) | (corr < {neg_cutoff})")
    else:
        raise ValueError(
            f'direction need to be in ["+", "-", "both"], got {direction}.'
        )

    print(f"{sum(selected)} correlation pairs selected.")

    corr_table = pd.DataFrame(
        {
            region_dim_a: pd.Series(total_results.obs_names)
            .iloc[total_results.X.row[selected]]
            .values,
            region_dim_b: pd.Series(total_results.var_names)
            .iloc[total_results.X.col[selected]]
            .values,
            "corr": total_results.X.data[selected],
        }
    )
    corr_table[f"{region_dim_a}_pos"] = corr_table[region_dim_a].map(
        total_results.obs["pos"]
    )
    corr_table[f"{region_dim_b}_pos"] = corr_table[region_dim_b].map(
        total_results.var["pos"]
    )
    return corr_table
