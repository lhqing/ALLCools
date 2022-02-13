from sklearn.decomposition import PCA
from dask_ml.decomposition import IncrementalPCA as dIPCA
import numpy as np


def simple_cca(adata, group_col, n_components=50, random_state=0):
    groups = adata.obs[group_col].unique()
    if len(groups) != 2:
        raise ValueError(
            f"CCA only handle 2 groups, adata.obs[{group_col}] has {len(groups)} different groups."
        )
    group_a, group_b = groups

    pca = PCA(n_components=n_components,
              copy=True,
              whiten=False,
              svd_solver='auto',
              tol=0.0,
              iterated_power='auto',
              random_state=random_state)

    a = adata[adata.obs[group_col] == group_a, :].X
    b = adata[adata.obs[group_col] == group_b, :].X
    X = a.dot(b.T)
    pc = pca.fit_transform(X)
    loading = pca.components_
    total_cc = np.concatenate([pc, loading.T], axis=0)
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
