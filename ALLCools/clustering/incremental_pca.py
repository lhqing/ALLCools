import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import IncrementalPCA as _IncrementalPCA
from ..count_matrix.zarr import dataset_to_array


def _normalize_per_cell(matrix, cell_sum):
    """Normalize matrix row sum to to cell_sum"""
    print("normalize per cell to CPM")
    if cell_sum is None:
        norm_vec = (matrix.sum(axis=1) + 1) / 1000000
    else:
        norm_vec = cell_sum / 1000000
        norm_vec = norm_vec.values
    norm_vec = norm_vec.astype(np.float32)
    matrix /= norm_vec[:, None]
    return matrix


class IncrementalPCA:
    def __init__(
        self,
        n_components=100,
        sparse=False,
        normalize_per_cell=True,
        log1p=True,
        scale=True,
        **kwargs,
    ):
        """
        Perform PCA for huge dataset that exceeds physical memory. Start from the raw count matrix stored in

        Parameters
        ----------
        n_components
            number of PCs to calculate
        sparse
            Whether treat the matrix as sparse matrix. If true, will
            1) load data as sparse matrix;
            2) do not perform mean center when scaling, only scale the std
        normalize_per_cell
            Normalize the matrix per cell
        log1p
            whether perform log1p transform
        scale
            whether scale the features
        kwargs
            parameters of sklearn.decomposition.IncrementalPCA
        """
        self.pca = _IncrementalPCA(n_components=n_components, **kwargs)
        self.sparse = sparse
        self.normalize_per_cell = normalize_per_cell
        self.log1p = log1p
        self.scale = scale
        self.scaler = None
        self.cell_sum = None
        self.use_features = None
        self.obs_dim = None
        self.var_dim = None
        self.load_chunk = None
        self._fit = False
        return

    def fit(
        self,
        ds,
        use_cells=None,
        use_features=None,
        chunk=500000,
        cell_sum=None,
        var_dim="gene",
        obs_dim="cell",
        load_chunk=None,
        random_shuffle=True,
    ):
        self.cell_sum = cell_sum
        self.use_features = use_features
        self.obs_dim = obs_dim
        self.var_dim = var_dim
        self.load_chunk = chunk if load_chunk is None else load_chunk

        # prepare index
        cell_index = ds.get_index(obs_dim)
        if use_cells is not None:
            cell_index = cell_index[cell_index.isin(use_cells)].copy()

        # random shuffle to make fitting more stable
        if random_shuffle:
            cell_order = cell_index.tolist()
            np.random.shuffle(cell_order)
            cell_order = pd.Index(cell_order)
        else:
            cell_order = cell_index

        # fit by chunks
        chunk_stds = []
        chunk_means = []
        for chunk_start in range(0, cell_order.size, chunk):
            print(f"Fitting {chunk_start}-{chunk_start + chunk}")
            _chunk_cells = cell_order[chunk_start : chunk_start + chunk]
            _chunk_matrix, _chunk_cells, _chunk_genes = dataset_to_array(
                ds,
                use_cells=_chunk_cells,
                use_genes=use_features,
                sparse=self.sparse,
                obs_dim=obs_dim,
                var_dim=var_dim,
                chunk=self.load_chunk,
            )
            if cell_sum is not None:
                _chunk_cell_sum = cell_sum.loc[_chunk_cells]
            else:
                _chunk_cell_sum = None
            _chunk_matrix = _chunk_matrix.astype(np.float32)

            # normalize cell counts
            if self.normalize_per_cell:
                _chunk_matrix = _normalize_per_cell(
                    matrix=_chunk_matrix, cell_sum=_chunk_cell_sum
                )

            # log transfer
            if self.log1p:
                print("log1p transform")
                _chunk_matrix = np.log1p(_chunk_matrix)

            # scale
            if self.scale:
                print("Scale")
                if self.scaler is None:
                    # assume the chunk is large enough, so only use the first chunk to fit
                    # e.g., 5,000,000 cells
                    self.scaler = StandardScaler(with_mean=not self.sparse)
                    _chunk_matrix = self.scaler.fit_transform(_chunk_matrix)
                else:
                    # transform remaining cells
                    _chunk_matrix = self.scaler.transform(_chunk_matrix)

            # save chunk stats for checking robustness
            chunk_stds.append(_chunk_matrix.std(axis=0))
            chunk_means.append(_chunk_matrix.mean(axis=0))

            # fit IncrementalPCA
            print("Fit PCA")
            self.pca.partial_fit(_chunk_matrix)

        self._fit = True
        return

    def transform(self, ds, use_cells=None, chunk=100000):
        if not self._fit:
            raise ValueError("fit first before transform")

        cell_index = ds.get_index(self.obs_dim)
        if use_cells is not None:
            cell_index = cell_index[cell_index.isin(use_cells)].copy()

        total_pcs = []
        for chunk_start in range(0, cell_index.size, chunk):
            print(f"Transforming {chunk_start}-{chunk_start + chunk}")
            _chunk_cells = cell_index[chunk_start : chunk_start + chunk]
            _chunk_matrix, _chunk_cells, _chunk_genes = dataset_to_array(
                ds,
                use_cells=_chunk_cells,
                use_genes=self.use_features,
                sparse=self.sparse,
                obs_dim=self.obs_dim,
                var_dim=self.var_dim,
                chunk=self.load_chunk,
            )
            if self.cell_sum is not None:
                _chunk_cell_sum = self.cell_sum.loc[_chunk_cells]
            else:
                _chunk_cell_sum = None
            _chunk_matrix = _chunk_matrix.astype(np.float32)

            # normalize cell counts
            if self.normalize_per_cell:
                _chunk_matrix = _normalize_per_cell(
                    matrix=_chunk_matrix, cell_sum=_chunk_cell_sum
                )

            # log transfer
            if self.log1p:
                print("log1p transform")
                _chunk_matrix = np.log1p(_chunk_matrix)

            # scale
            if self.scale:
                print("Scale")
                if self.scaler is None:
                    # this shouldn't happen in transform
                    raise ValueError("scale is True, but scaler not exist")
                else:
                    # transform remaining cells
                    _chunk_matrix = self.scaler.transform(_chunk_matrix)

            # transform
            print("Transform PCA")
            pcs = self.pca.transform(_chunk_matrix)
            pcs = pd.DataFrame(pcs, index=_chunk_cells)
            total_pcs.append(pcs)
        total_pcs = pd.concat(total_pcs)
        return total_pcs

    def fit_transform(
        self,
        ds,
        use_cells=None,
        use_features=None,
        chunk=500000,
        cell_sum=None,
        var_dim="gene",
        obs_dim="cell",
        load_chunk=None,
        random_shuffle=True,
    ):
        self.fit(
            ds,
            use_cells=use_cells,
            use_features=use_features,
            chunk=chunk,
            cell_sum=cell_sum,
            var_dim=var_dim,
            obs_dim=obs_dim,
            load_chunk=load_chunk,
            random_shuffle=random_shuffle,
        )

        total_pcs = self.transform(ds, use_cells=use_cells, chunk=self.load_chunk)
        return total_pcs
