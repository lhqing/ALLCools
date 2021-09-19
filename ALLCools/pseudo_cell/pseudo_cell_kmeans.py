import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from scipy.sparse import issparse, csr_matrix, vstack
import anndata


def _kmeans_iter(matrix, cells, max_group_size, k=None, parent_label=None):
    """Iteratively split matrix using kmeans until all the group size < max_group_size"""
    if max_group_size <= 1:
        return pd.Series(range(cells.size), index=cells).astype(str)
    if k is None:
        k = matrix.shape[0] // max_group_size + 1
    if matrix.shape[0] <= k:
        return pd.Series(np.ones_like(cells), index=cells).astype(str)

    # iteratively split cells until all the
    mbk = MiniBatchKMeans(n_clusters=k,
                          init='k-means++',
                          max_iter=100,
                          batch_size=100,
                          verbose=0,
                          compute_labels=True,
                          random_state=0,
                          tol=0.0,
                          max_no_improvement=10,
                          init_size=3 * k,
                          n_init=5,
                          reassignment_ratio=0.1)
    mbk.fit(matrix)
    labels = pd.Series(mbk.labels_, index=cells).astype(str)
    if parent_label is not None:
        labels = labels.apply(lambda i: f'{parent_label}|{i}')

    records = []
    for label, group in labels.groupby(labels):
        if group.size <= max_group_size:
            records.append(group)
        else:
            # group is over sized, keep splitting
            # get matrix of this group
            group_matrix = matrix[cells.isin(group.index), :]
            group_cells = group.index
            new_k = group_matrix.shape[0] // max_group_size + 1
            records += _kmeans_iter(matrix=group_matrix,
                                    cells=group_cells,
                                    k=new_k,
                                    max_group_size=max_group_size,
                                    parent_label=label)
    if parent_label is None:
        # the top function call
        records = pd.concat(records)
        rename_cluster = {c: str(i) for i, c in enumerate(records.value_counts().index)}
        records = records.map(rename_cluster).sort_index()
    return records


def _calculate_cell_group(clusters, total_matrix, cluster_size_cutoff=100, max_group_size=25):
    cluster_cells = clusters.value_counts()
    max_group_sizes = (cluster_cells // cluster_size_cutoff + 1).astype(int)
    max_group_sizes[max_group_sizes > max_group_size] = max_group_size

    records = []
    for cluster, max_group_size in max_group_sizes.items():
        cells = clusters[clusters == cluster].index
        matrix = total_matrix[clusters == cluster].copy()
        record = _kmeans_iter(matrix, cells, max_group_size, k=None, parent_label=None)
        record = record.apply(lambda i: f'{cluster}::{i}')
        records.append(record)
    records = pd.concat(records)
    return records


def _merge_pseudo_cell(adata, aggregate_func='downsample'):
    is_sparse = issparse(adata.X)
    # merge ia to balanced ia (bia)
    balanced_matrix = []
    obs = {}
    for group, cells in adata.obs_names.groupby(adata.obs['cell_group']).items():
        n_cells = cells.size
        if cells.size == 1:
            group_data = adata.var_vector(cells[0]).ravel()
        else:
            if aggregate_func == 'sum':
                try:
                    group_data = adata[cells, :].X.sum(axis=0).A1
                except AttributeError:
                    group_data = np.array(adata[cells, :].X.sum(axis=0)).ravel()
            elif aggregate_func == 'mean':
                try:
                    group_data = adata[cells, :].X.mean(axis=0).A1
                except AttributeError:
                    group_data = np.array(adata[cells, :].X.mean(axis=0)).ravel()
            elif aggregate_func == 'median':
                try:
                    group_data = adata[cells, :].X.median(axis=0).A1
                except AttributeError:
                    group_data = np.array(adata[cells, :].X.median(axis=0)).ravel()
            elif aggregate_func == 'downsample':
                cell = np.random.choice(cells, 1)
                try:
                    group_data = adata[cell, :].X.sum(axis=0).A1  # just to reduce dim
                except AttributeError:
                    group_data = np.array(adata[cell, :].X).ravel()
            else:
                raise ValueError(f'aggregate_func can only be ["sum", "mean", "median"], got "{aggregate_func}"')
        balanced_matrix.append(csr_matrix(group_data))
        obs[group] = {'n_cells': n_cells}
    balanced_matrix = vstack(balanced_matrix)
    if not is_sparse:
        balanced_matrix = balanced_matrix.toarray()

    pseudo_cell_adata = anndata.AnnData(balanced_matrix,
                                        obs=pd.DataFrame(obs).T,
                                        var=adata.var.copy())
    return pseudo_cell_adata


def generate_pseudo_cells(adata,
                          cluster_col='leiden',
                          obsm='X_pca',
                          cluster_size_cutoff=100,
                          max_group_size=25,
                          aggregate_func='downsample'):
    # determine cell group
    clusters = adata.obs[cluster_col]
    total_matrix = adata.obsm[obsm]
    cell_group = _calculate_cell_group(clusters=clusters,
                                       total_matrix=total_matrix,
                                       cluster_size_cutoff=cluster_size_cutoff,
                                       max_group_size=max_group_size)
    adata.obs['cell_group'] = cell_group
    pseudo_cell_adata = _merge_pseudo_cell(adata=adata,
                                           aggregate_func=aggregate_func)
    pseudo_cell_adata.obs[cluster_col] = pseudo_cell_adata.obs_names.str.split('::').str[0]
    return pseudo_cell_adata
