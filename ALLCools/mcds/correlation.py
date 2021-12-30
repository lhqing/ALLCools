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
            out[idx] = (a[i] - mu_a[i]) @ (b[j] -
                                           mu_b[j]) / k / _sig_a / _sig_b
    return out


def _corr_preprocess(da, region_dim, sample_mch, sample_mcg, cpu=1):
    df = da.to_pandas()
    adata = anndata.AnnData(X=df.values,
                            obs=pd.DataFrame([], index=df.index),
                            var=pd.DataFrame([], index=df.columns))

    start = da.coords[f'{region_dim}_start']
    end = da.coords[f'{region_dim}_end']
    center = (end + start) // 2
    adata.var['center'] = center

    # standardize global
    sample_mcg = (sample_mcg - sample_mcg.mean()) / sample_mcg.std()
    sample_mch = (sample_mch - sample_mch.mean()) / sample_mch.std()
    # add global
    adata.obs['sample_mch'] = sample_mch
    adata.obs['sample_mcg'] = sample_mcg

    # remove dmr that has 0 std
    # otherwise regress out will fail.
    adata = adata[:, adata.X.std(axis=0) > 0].copy()

    # regress out global mCH and mCG
    # happens on each regions
    sc.pp.regress_out(adata, keys=['sample_mch', 'sample_mcg'], n_jobs=cpu)
    sc.pp.scale(adata)
    return adata


def corr(adata_a, adata_b, max_dist, cpu=1, method='pearson', null=None, null_n=1e6):
    # create mask
    mask = coo_matrix(
        np.abs((adata_a.var['center'].values[:, None] -
                adata_b.var['center'].values[None, :])) < max_dist)

    if null == 'sample':
        obs_names = adata_a.obs_names.values.copy()
        np.random.shuffle(obs_names)
        adata_a.obs_names = obs_names
        adata_a = adata_a[adata_b.obs_names, :].copy()  # reorder adata_a so its shuffled

        if mask.data.size > null_n:
            rand_loc = np.random.choice(range(mask.data.size), size=int(null_n), replace=False)
            rand_loc = sorted(rand_loc)
            mask = coo_matrix((mask.data[rand_loc],
                               (mask.row[rand_loc], mask.col[rand_loc])),
                              shape=mask.shape)
    elif null == 'region':
        row, col = mask.shape
        null_ratio = min(1, null_n / row / col)

        null_mask = np.random.rand(*mask.shape) < null_ratio
        mask = coo_matrix(null_mask & ~mask.toarray())
    else:
        pass

    # adata is sample by region
    # _corr input is region by sample
    a = adata_a.X.T
    b = adata_b.X.T

    if method.lower()[0] == 'p':
        pass
    elif method.lower()[0] == 's':
        # turn a, b in to rank matrix
        a = a.argsort(axis=1).argsort(axis=1)
        b = b.argsort(axis=1).argsort(axis=1)
    else:
        raise ValueError('Method can only be pearson or spearman')

    a = a.astype(np.float32)
    b = b.astype(np.float32)

    chunk_size = mask.data.size // cpu + 1

    # trigger njit
    _corr(np.array([[0, 1]]), np.array([[0, 1]]),
          np.array([0]), np.array([0]))

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for i, chunk_start in enumerate(range(0, mask.data.size, chunk_size)):
            future = exe.submit(
                _corr,
                a=a,
                b=b,
                row=mask.row[chunk_start:chunk_start + chunk_size],
                col=mask.col[chunk_start:chunk_start + chunk_size])
            futures[future] = i

        total_data = {}
        for future in as_completed(futures):
            chunk_id = futures[future]
            total_data[chunk_id] = future.result()
        mask.data = np.concatenate([total_data[k] for k in sorted(total_data.keys())])

    result = anndata.AnnData(X=mask.tocsr(),
                             obs=pd.DataFrame([], index=adata_a.var_names),
                             var=pd.DataFrame([], index=adata_b.var_names))
    return result


def region_correlation(data_a, data_b, sample_mch, sample_mcg, method='pearson', max_dist=10000000,
                       cpu=1, null=None, null_n=1e6):
    region_dim_a = data_a.dims[1]
    region_dim_b = data_b.dims[1]

    total_results = []
    for chrom_a, chrom_da_a in data_a.groupby(f"{region_dim_a}_chrom"):
        for chrom_b, chrom_da_b in data_b.groupby(f"{region_dim_b}_chrom"):
            if chrom_a != chrom_b:
                continue
            print(f'Calculating {chrom_a}')
            adata_a = _corr_preprocess(da=chrom_da_a,
                                       region_dim=region_dim_a,
                                       sample_mch=sample_mch,
                                       sample_mcg=sample_mcg,
                                       cpu=cpu)
            adata_b = _corr_preprocess(da=chrom_da_b,
                                       region_dim=region_dim_b,
                                       sample_mch=sample_mch,
                                       sample_mcg=sample_mcg,
                                       cpu=cpu)
            corr_matrix = corr(adata_a,
                               adata_b,
                               max_dist=max_dist,
                               method=method,
                               null=null,
                               null_n=null_n,
                               cpu=cpu)
            total_results.append(corr_matrix)

    total_results = anndata.concat(total_results)
    return total_results
