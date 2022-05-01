import numpy as np
from dask_ml.decomposition import IncrementalPCA as dIPCA
from scipy.stats import zscore
from sklearn.decomposition import PCA


def cca(data1, data2, scale1=False, scale2=False, n_components=50, random_state=42):
    if scale1:
        data1 = zscore(data1, axis=1)
    if scale2:
        data2 = zscore(data2, axis=1)

    X = data1.dot(data2.T)
    pca = PCA(n_components=n_components,
              copy=True,
              whiten=False,
              svd_solver='auto',
              tol=0.0,
              iterated_power='auto',
              random_state=random_state)
    pc = pca.fit_transform(X)
    loading = pca.components_
    return pc, loading.T


def adata_cca(adata, group_col, separate_scale=True, n_components=50, random_state=42):
    groups = adata.obs[group_col].unique()
    if len(groups) != 2:
        raise ValueError(
            f"CCA only handle 2 groups, adata.obs[{group_col}] has {len(groups)} different groups."
        )
    group_a, group_b = groups
    a = adata[adata.obs[group_col] == group_a, :].X
    b = adata[adata.obs[group_col] == group_b, :].X

    pc, loading = cca(data1=a,
                      data2=b,
                      scale1=separate_scale,
                      scale2=separate_scale,
                      n_components=n_components,
                      random_state=random_state)
    total_cc = np.concatenate([pc, loading], axis=0)
    adata.obsm['X_cca'] = total_cc
    return


def incremental_cca(a, b, max_chunk_size=10000, random_state=0):
    """
    Perform Incremental CCA by chunk dot product and IncrementalPCA

    Parameters
    ----------
    a
        dask.Array of dataset a
    b
        dask.Array of dataset b
    max_chunk_size
        Chunk size for Incremental fit and transform, the larger the better as long as MEM is enough
    random_state

    Returns
    -------
    Top CCA components
    """
    pca = dIPCA(n_components=50,
                whiten=False,
                copy=True,
                batch_size=None,
                svd_solver='auto',
                iterated_power=0,
                random_state=random_state)

    # partial fit
    n_sample = a.shape[0]
    n_chunks = n_sample // max_chunk_size + 1
    chunk_size = int(n_sample / n_chunks) + 1
    for chunk_start in range(0, n_sample, chunk_size):
        print(chunk_start)
        X_chunk = a[chunk_start:chunk_start + chunk_size, :].dot(b.T)
        pca.partial_fit(X_chunk)

    # transform
    pcs = []
    for chunk_start in range(0, n_sample, chunk_size):
        print(chunk_start)
        X_chunk = a[chunk_start:chunk_start + chunk_size, :].dot(b.T)
        pc_chunk = pca.transform(X_chunk).compute()
        pcs.append(pc_chunk)
    pcs = np.concatenate(pcs)

    # concatenate CCA
    total_cc = np.concatenate([pcs, pca.components_.T])
    return total_cc


def _tf_idf(data, col_sum, scale_factor):
    row_sum = data.sum(axis=1).A1.astype(int)
    data.data = data.data / np.repeat(row_sum, row_sum)
    data.data = np.log(data.data * scale_factor + 1)
    idf = np.log(1 + data.shape[0] / col_sum)
    return idf


def lsi_cca(data1, data2, min_cov=5, scale_factor=100000, n_components=50, random_state=42):
    col_sum1 = data1.sum(axis=0).A1
    col_sum2 = data2.sum(axis=0).A1
    binfilter = np.logical_and(col_sum1 > min_cov, col_sum2 > min_cov)
    data1 = data1[:, binfilter]
    data2 = data2[:, binfilter]

    idf1 = _tf_idf(data1, col_sum1, scale_factor)
    idf2 = _tf_idf(data2, col_sum2, scale_factor)
    tf = data1.multiply(idf1).dot(data2.multiply(idf2).T)

    pca = PCA(n_components=n_components,
              copy=True,
              whiten=False,
              svd_solver='auto',
              tol=0.0,
              iterated_power='auto',
              random_state=random_state)
    pc = pca.fit_transform(tf)
    loading = pca.components_
    return pc, loading.T
