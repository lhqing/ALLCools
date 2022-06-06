import numpy as np
from dask_ml.decomposition import IncrementalPCA as dIPCA
from scipy.stats import zscore
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.extmath import safe_sparse_dot
from ..clustering.lsi import tf_idf


def cca(data1, data2, scale1=True, scale2=True, n_components=50, max_cc_cell=20000, random_state=0):
    if scale1:
        data1 = zscore(data1, axis=1)
    if scale2:
        data2 = zscore(data2, axis=1)
    np.random.seed(random_state)
    model = TruncatedSVD(n_components=n_components, algorithm='arpack', random_state=random_state)

    U = model.fit_transform(data1.dot(data2.T))
    sel_dim = (model.singular_values_ != 0)
    if max_cc_cell > data2.shape[0]:
        V = model.components_[sel_dim].T
    else:
        V = safe_sparse_dot(safe_sparse_dot(U.T[sel_dim], data1), data2.T).T
        V = V / np.square(model.singular_values_)
    if max_cc_cell > data1.shape[0]:
        U = U[:, sel_dim] / model.singular_values_
    else:
        U = safe_sparse_dot(data1, safe_sparse_dot(model.components_[sel_dim], data2).T)
        U = U / model.singular_values_
    return U, V


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
    raise NotImplementedError
    # TODO PC is wrong
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


def lsi_cca(data1, data2, scale_factor=100000, n_components=50, max_cc_cell=20000, chunk_size=50000):
    np.random.seed(0)
    sel1 = np.random.choice(np.arange(data1.shape[0]), min(max_cc_cell, data1.shape[0]), False)
    sel2 = np.random.choice(np.arange(data2.shape[0]), min(max_cc_cell, data2.shape[0]), False)
    tf1, idf1 = tf_idf(data1[sel1], scale_factor=scale_factor)
    tf2, idf2 = tf_idf(data2[sel2], scale_factor=scale_factor)
    tf = tf1.dot(tf2.T)
    model = TruncatedSVD(n_components=n_components, algorithm='arpack', random_state=0)
    U = model.fit_transform(tf)
    seldim = (model.singular_values_ != 0)
    if max_cc_cell > data2.shape[0]:
        V = model.components_[seldim].T
    else:
        V = np.concatenate([safe_sparse_dot(safe_sparse_dot(U.T[seldim], tf1),
                                            tf_idf(data2[chunk_start:(chunk_start + chunk_size)],
                                                   scale_factor=scale_factor, idf=idf2)[0].T).T
                            for chunk_start in np.arange(0, data2.shape[0], chunk_size)], axis=0)
        V = V / np.square(model.singular_values_)
    if max_cc_cell > data1.shape[0]:
        U = U[:, seldim] / model.singular_values_
    else:
        U = np.concatenate([safe_sparse_dot(
            tf_idf(data1[chunk_start:(chunk_start + chunk_size)], scale_factor=scale_factor, idf=idf1)[0],
            safe_sparse_dot(model.components_[seldim], tf2).T)
            for chunk_start in np.arange(0, data1.shape[0], chunk_size)], axis=0)
        U = U / model.singular_values_
    return U, V
