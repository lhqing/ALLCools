import dask
import numpy as np
import pandas as pd
import anndata
import xarray as xr
from scipy.sparse import csr_matrix, vstack


def open_zarr(path, obs_dim="cell"):
    ds = xr.open_mfdataset(path, engine="zarr", concat_dim=obs_dim, combine="nested")
    return ds


def dataset_to_array(
    ds,
    use_cells=None,
    use_genes=None,
    sparse=True,
    obs_dim="cell",
    var_dim="gene",
    chunk=100000,
):
    cell_index = ds.get_index(obs_dim)
    gene_index = ds.get_index(var_dim)

    if use_cells is None:
        use_cells = cell_index
    else:
        use_cells = cell_index[cell_index.isin(use_cells)]
    if use_genes is None:
        use_genes = gene_index
    else:
        use_genes = gene_index[gene_index.isin(use_genes)]

    # load data by chunk
    data = []
    for chunk_start in range(0, use_cells.size, chunk):
        _chunk_cells = use_cells[chunk_start : chunk_start + chunk]
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            _chunk_ds = ds[f"{var_dim}_da"].sel(
                {obs_dim: _chunk_cells, var_dim: use_genes}
            )
        if sparse:
            data.append(csr_matrix(_chunk_ds.values))
        else:
            data.append(_chunk_ds.values)
    if sparse:
        data = vstack(data)
    else:
        data = np.concatenate(data)
    return data, use_cells, use_genes


def dataset_to_adata(
    ds,
    use_cells=None,
    use_genes=None,
    sparse=True,
    obs_dim="cell",
    var_dim="gene",
    chunk=100000,
):
    data, use_cells, use_genes = dataset_to_array(
        ds,
        use_cells=use_cells,
        use_genes=use_genes,
        sparse=sparse,
        obs_dim=obs_dim,
        var_dim=var_dim,
        chunk=chunk,
    )

    adata = anndata.AnnData(
        X=data,
        obs=pd.DataFrame([], index=use_cells),
        var=pd.DataFrame([], index=use_genes),
    )
    return adata
