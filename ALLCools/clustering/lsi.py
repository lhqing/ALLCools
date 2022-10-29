import joblib
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.validation import check_is_fitted


def tf_idf(data, scale_factor=100000, idf=None):
    sparse_input = issparse(data)

    if idf is None:
        # add small value in case down sample creates empty feature
        _col_sum = data.sum(axis=0)
        if sparse_input:
            col_sum = _col_sum.A1.astype(np.float32) + 0.00001
        else:
            col_sum = _col_sum.ravel().astype(np.float32) + 0.00001
        idf = np.log(1 + data.shape[0] / col_sum).astype(np.float32)
    else:
        idf = idf.astype(np.float32)

    _row_sum = data.sum(axis=1)
    if sparse_input:
        row_sum = _row_sum.A1.astype(np.float32) + 0.00001
    else:
        row_sum = _row_sum.ravel().astype(np.float32) + 0.00001

    tf = data.astype(np.float32)

    if sparse_input:
        tf.data = tf.data / np.repeat(row_sum, tf.getnnz(axis=1))
        tf.data = np.log1p(np.multiply(tf.data, scale_factor, dtype="float32"))
        tf = tf.multiply(idf)
    else:
        tf = tf / row_sum[:, np.newaxis]
        tf = np.log1p(np.multiply(tf, scale_factor, dtype="float32"))
        tf = tf * idf
    return tf, idf


def lsi(
    adata,
    scale_factor=100000,
    n_components=100,
    algorithm="arpack",
    obsm="X_pca",
    random_state=0,
    fit_size=None,
):
    """
    Run TF-IDF on the binarized adata.X, followed by TruncatedSVD and then scale the components by singular values.

    Parameters
    ----------
    adata
        AnnData object
    scale_factor
        scale factor for TF-IDF
    n_components
        number of components to keep
    algorithm
        algorithm to use for TruncatedSVD
    obsm
        key in adata.obsm to store the components in
    random_state
        random state for reproducibility
    fit_size
        Ratio or absolute int value, use to downsample when fitting the SVD to speed up run time.
    """
    # Raise deprecationWarning
    from warnings import warn

    warn("lsi function is deprecated, use LSI class instead", DeprecationWarning)

    # tf-idf
    data = adata.X.astype(np.int8).copy()
    tf, _ = tf_idf(data, scale_factor)
    n_rows, n_cols = tf.shape
    n_components = min(n_rows, n_cols, n_components)
    svd = TruncatedSVD(n_components=n_components, algorithm=algorithm, random_state=random_state)

    if fit_size is None:
        # fit the SVD using all rows
        matrix_reduce = svd.fit_transform(tf)
    elif fit_size >= n_rows:
        # fit size is larger than actual data size
        matrix_reduce = svd.fit_transform(tf)
    else:
        # fit the SVD using partial rows to speed up
        if fit_size < 1:
            fit_size = max(int(n_rows * fit_size), n_components)
        use_cells = pd.Series(range(n_rows)).sample(fit_size, random_state=random_state).sort_index().tolist()
        svd.fit(tf.tocsr()[use_cells, :])
        matrix_reduce = svd.transform(tf)

    matrix_reduce = matrix_reduce / svd.singular_values_

    # PCA is the default name for many following steps in scanpy, use the name here for convenience.
    # However, this is not PCA
    adata.obsm[obsm] = matrix_reduce
    return svd


class LSI:
    def __init__(
        self,
        scale_factor=100000,
        n_components=100,
        algorithm="arpack",
        random_state=0,
        idf=None,
        model=None,
    ):
        self.scale_factor = scale_factor
        if idf is not None:
            self.idf = idf.copy()
        if model is not None:
            self.model = model
        else:
            self.model = TruncatedSVD(n_components=n_components, algorithm=algorithm, random_state=random_state)
        self.random_state = random_state

    def _downsample_data(self, data, downsample):
        np.random.seed(self.random_state)
        if downsample is not None and downsample < data.shape[0]:
            use_row_idx = np.sort(np.random.choice(np.arange(0, data.shape[0]), downsample, replace=False))
            data = data[use_row_idx, :]
        return data

    @staticmethod
    def _get_data(data):
        if isinstance(data, AnnData):
            data = data.X
        return data

    def fit(self, data, downsample=None):
        data = self._get_data(data)
        data = self._downsample_data(data, downsample)
        tf, idf = tf_idf(data, self.scale_factor)
        self.idf = idf.copy()
        n_rows, n_cols = tf.shape
        self.model.n_components = min(n_rows - 1, n_cols - 1, self.model.n_components)
        self.model.fit(tf)
        return self

    def fit_transform(self, data, downsample=None, obsm_name="X_lsi"):
        _data = self._get_data(data)
        _data = self._downsample_data(_data, downsample)
        tf, idf = tf_idf(_data, self.scale_factor)
        self.idf = idf.copy()
        n_rows, n_cols = tf.shape
        self.model.n_components = min(n_rows - 1, n_cols - 1, self.model.n_components)
        tf_reduce = self.model.fit_transform(tf)

        tf_reduce = tf_reduce / self.model.singular_values_

        if isinstance(data, AnnData):
            data.obsm[obsm_name] = tf_reduce
        else:
            return tf_reduce

    def transform(self, data, chunk_size=50000, obsm_name="X_lsi"):
        _data = self._get_data(data)

        check_is_fitted(self.model)
        tf_reduce = []
        for chunk_start in np.arange(0, _data.shape[0], chunk_size):
            tf, _ = tf_idf(_data[chunk_start : (chunk_start + chunk_size)], self.scale_factor, self.idf)
            tf_reduce.append(self.model.transform(tf))

        tf_reduce = np.concatenate(tf_reduce, axis=0)

        tf_reduce = tf_reduce / self.model.singular_values_

        if isinstance(data, AnnData):
            data.obsm[obsm_name] = tf_reduce
        else:
            return tf_reduce

    def save(self, path):
        joblib.dump(self, path)
