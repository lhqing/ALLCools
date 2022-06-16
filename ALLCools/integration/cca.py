import numpy as np
# from dask_ml.decomposition import IncrementalPCA as dIPCA
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.extmath import safe_sparse_dot
from ..clustering.lsi import tf_idf
from sklearn.preprocessing import StandardScaler
from scipy.sparse import issparse


def cca(data1,
        data2,
        scale1=True,
        scale2=True,
        n_components=50,
        max_cc_cell=20000,
        chunk_size=50000,
        random_state=0):
    np.random.seed(random_state)

    # downsample cells
    if max_cc_cell < data1.shape[0]:
        sel1 = np.sort(
            np.random.choice(
                np.arange(data1.shape[0]),
                min(max_cc_cell, data1.shape[0]),
                False
            )
        )
        tf_data1 = data1[sel1, :]
    else:
        tf_data1 = data1
    if max_cc_cell < data2.shape[0]:
        sel2 = np.sort(
            np.random.choice(
                np.arange(data2.shape[0]),
                min(max_cc_cell, data2.shape[0]),
                False
            )
        )
        tf_data2 = data2[sel2, :]
    else:
        tf_data2 = data2

    # sparse matrix to dense
    if issparse(tf_data1):
        tf_data1 = tf_data1.toarray()
    if issparse(tf_data2):
        tf_data2 = tf_data2.toarray()

    # scale input features
    scaler1 = StandardScaler()
    if scale1:
        tf_data1 = scaler1.fit_transform(tf_data1)
    scaler2 = StandardScaler()
    if scale2:
        tf_data2 = scaler2.fit_transform(tf_data2)

    # CCA decomposition
    model = TruncatedSVD(n_components=n_components,
                         algorithm='arpack',
                         random_state=random_state)
    U = model.fit_transform(tf_data1.dot(tf_data2.T))

    # select dimensions with non-zero singular values
    sel_dim = model.singular_values_ != 0
    print('non zero dims', sel_dim.sum())
    nnz_components = model.components_[:, sel_dim]
    nnz_singular_values = model.singular_values_[sel_dim]

    if max_cc_cell > data2.shape[0]:
        V = nnz_components.T
    else:
        V = []
        for chunk_start in np.arange(0, data2.shape[0], chunk_size):
            if issparse(data2):
                tmp = data2[chunk_start:(chunk_start + chunk_size)].toarray()
            else:
                tmp = data2[chunk_start:(chunk_start + chunk_size)]
            if scale2:
                tmp = scaler2.transform(tmp)
            # calculate V by fitted U and the sampled data1
            V.append(np.dot(np.dot(U.T[sel_dim], tf_data1), tmp.T).T)
        V = np.concatenate(V, axis=0)
        V = V / np.square(nnz_singular_values)

    if max_cc_cell > data1.shape[0]:
        U = U[:, sel_dim] / nnz_singular_values
    else:
        U = []
        for chunk_start in np.arange(0, data1.shape[0], chunk_size):
            if issparse(data1):
                tmp = data1[chunk_start:(chunk_start + chunk_size)].toarray()
            else:
                tmp = data1[chunk_start:(chunk_start + chunk_size)]
            if scale1:
                tmp = scaler1.transform(tmp)
            U.append(np.dot(tmp, np.dot(nnz_components, tf_data2).T))
        U = np.concatenate(U, axis=0)
        U = U / nnz_singular_values
    return U, V


def adata_cca(adata, group_col, separate_scale=True, n_components=50, random_state=42):
    groups = adata.obs[group_col].unique()
    if len(groups) != 2:
        raise ValueError(
            f"CCA only handle 2 groups, "
            f"adata.obs[{group_col}] has {len(groups)} different groups."
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


# def incremental_cca(a, b, max_chunk_size=10000, random_state=0):
#     """
#     Perform Incremental CCA by chunk dot product and IncrementalPCA
#
#     Parameters
#     ----------
#     a
#         dask.Array of dataset a
#     b
#         dask.Array of dataset b
#     max_chunk_size
#         Chunk size for Incremental fit and transform, the larger the better as long as MEM is enough
#     random_state
#
#     Returns
#     -------
#     Top CCA components
#     """
#     raise NotImplementedError
#     # TODO PC is wrong
#     pca = dIPCA(n_components=50,
#                 whiten=False,
#                 copy=True,
#                 batch_size=None,
#                 svd_solver='auto',
#                 iterated_power=0,
#                 random_state=random_state)
#
#     # partial fit
#     n_sample = a.shape[0]
#     n_chunks = n_sample // max_chunk_size + 1
#     chunk_size = int(n_sample / n_chunks) + 1
#     for chunk_start in range(0, n_sample, chunk_size):
#         print(chunk_start)
#         X_chunk = a[chunk_start:chunk_start + chunk_size, :].dot(b.T)
#         pca.partial_fit(X_chunk)
#
#     # transform
#     pcs = []
#     for chunk_start in range(0, n_sample, chunk_size):
#         print(chunk_start)
#         X_chunk = a[chunk_start:chunk_start + chunk_size, :].dot(b.T)
#         pc_chunk = pca.transform(X_chunk).compute()
#         pcs.append(pc_chunk)
#     pcs = np.concatenate(pcs)
#
#     # concatenate CCA
#     total_cc = np.concatenate([pcs, pca.components_.T])
#     return total_cc


def lsi_cca(data1,
            data2,
            scale_factor=100000,
            n_components=50,
            max_cc_cell=20000,
            chunk_size=50000,
            min_cov_filter=5):
    np.random.seed(0)

    # downsample data1 and data2 to run tf_idf and CCA
    if max_cc_cell < data1.shape[0]:
        sel1 = np.sort(
            np.random.choice(
                np.arange(data1.shape[0]),
                max_cc_cell,
                False
            )
        )
        tf_data1 = data1[sel1, :]
    else:
        tf_data1 = data1
    if max_cc_cell < data2.shape[0]:
        sel2 = np.sort(
            np.random.choice(
                np.arange(data2.shape[0]),
                max_cc_cell,
                False
            )
        )
        tf_data2 = data2[sel2, :]
    else:
        tf_data2 = data2

    # filter bin to make sure the min_cov_filter is satisfied
    col_sum1 = tf_data1.sum(axis=0).A1
    col_sum2 = tf_data2.sum(axis=0).A1
    # the same bin_filter will also be used
    # in the chunk transfer below
    bin_filter = np.logical_and(col_sum1 > min_cov_filter,
                                col_sum2 > min_cov_filter)
    tf1, idf1 = tf_idf(tf_data1[:, bin_filter], scale_factor=scale_factor)
    tf2, idf2 = tf_idf(tf_data2[:, bin_filter], scale_factor=scale_factor)

    # CCA part
    model = TruncatedSVD(n_components=n_components,
                         algorithm='arpack',
                         random_state=0)
    tf = tf1.dot(tf2.T)
    U = model.fit_transform(tf)

    # select non-zero singular values
    # transform the whole dataset
    sel_dim = (model.singular_values_ != 0)
    nnz_singular_values = model.singular_values_[sel_dim]
    nnz_components = model.components_[sel_dim]
    if max_cc_cell > data2.shape[0]:
        V = nnz_components.T
    else:
        # use the safe_sparse_dot to avoid memory error
        V = np.concatenate(
            [safe_sparse_dot(
                safe_sparse_dot(U.T[sel_dim], tf1),
                tf_idf(
                    data2[chunk_start:(chunk_start + chunk_size)][:, bin_filter],
                    scale_factor=scale_factor,
                    idf=idf2
                )[0].T  # [0] is the tf
            ).T
             for chunk_start in np.arange(0, data2.shape[0], chunk_size)],
            axis=0
        )
        V = V / np.square(nnz_singular_values)
    if max_cc_cell > data1.shape[0]:
        U = U[:, sel_dim] / nnz_singular_values
    else:
        U = np.concatenate(
            [safe_sparse_dot(
                tf_idf(
                    data1[chunk_start:(chunk_start + chunk_size)][:, bin_filter],
                    scale_factor=scale_factor,
                    idf=idf1
                )[0],  # [0] is the tf
                safe_sparse_dot(nnz_components, tf2).T
            )
                for chunk_start in np.arange(0, data1.shape[0], chunk_size)],
            axis=0
        )
        U = U / nnz_singular_values
    return U, V
