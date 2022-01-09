import h5py
import anndata
import pandas as pd
import numpy as np
import xarray as xr
import scipy.sparse as ss

SNAP_BD_FIELDS = {
    "TN": "TotalFrag",
    "UM": "UniqueMapFrag",
    "SE": "SingleEndFrag",
    "SA": "SecondaryAlign",
    "PE": "PairEndFrag",
    "PP": "ProperlyPairedFrag",
    "PL": "ProperFragLength",
    "US": "UsableFrag",
    "UQ": "UniqFrag",
    "CM": "chrMFrag",
}


def read_snap(file_path, bin_kind=5000):
    """Read snap hdf5 file into anndata.AnnData."""
    if bin_kind == "gene":
        return _read_snap_genes(file_path)
    elif isinstance(bin_kind, int):
        return _read_snap_bins(file_path, bin_size=bin_kind)
    else:
        raise NotImplementedError


def _read_snap_meta(f):
    cell_barcode = np.array(f["/BD/name"]).astype(str)
    cell_meta = pd.DataFrame([], index=cell_barcode)
    for name in f["/BD/"].keys():
        if name == "name":
            continue
        if name in SNAP_BD_FIELDS:
            col_name = SNAP_BD_FIELDS[name]
        else:
            col_name = name
        cell_meta[col_name] = np.array(f[f"/BD/{name}"])
    return cell_meta, cell_barcode


def _read_snap_genes(file_path):
    """Read gene data from snap hdf5 file into anndata.AnnData"""

    with h5py.File(file_path) as f:
        data = np.array(f[f"/GM/count"], dtype=np.uint8)
        idx = np.array(f[f"/GM/idx"]) - 1  # R is 1 based, python is 0 based
        idy = np.array(f[f"/GM/idy"]) - 1  # R is 1 based, python is 0 based
        gene_id = np.array(f[f"/GM/name"]).astype(str)

        cell_meta, cell_barcode = _read_snap_meta(f)
        data = ss.coo_matrix(
            (data, (idx, idy)), shape=(cell_barcode.size, gene_id.size), dtype=np.uint8
        ).tocsc()

    adata = anndata.AnnData(
        X=data, obs=cell_meta, var=pd.DataFrame(index=pd.Index(gene_id, name="gene"))
    )
    return adata


def _read_snap_bins(file_path, bin_size=5000):
    """Read bin data from snap hdf5 file into anndata.AnnData"""
    with h5py.File(file_path) as f:
        data = np.array(f[f"/AM/{bin_size}/count"], dtype=np.uint8)
        idx = np.array(f[f"/AM/{bin_size}/idx"]) - 1  # R is 1 based, python is 0 based
        idy = np.array(f[f"/AM/{bin_size}/idy"]) - 1  # R is 1 based, python is 0 based

        bin_chrom = np.array(f[f"/AM/{bin_size}/binChrom"]).astype(str)
        bin_start = np.array(f[f"/AM/{bin_size}/binStart"]) - 1  # 0 based bed format
        bin_end = (
            np.array(f[f"/AM/{bin_size}/binStart"]) - 1 + bin_size
        )  # 0 based bed format
        bin_id = np.core.defchararray.add(
            np.core.defchararray.add(bin_chrom, "-"), bin_start.astype(str)
        )

        cell_meta, cell_barcode = _read_snap_meta(f)
        data = ss.coo_matrix(
            (data, (idx, idy)), shape=(cell_barcode.size, bin_id.size), dtype=np.uint8
        ).tocsc()

    adata = anndata.AnnData(
        X=data,
        obs=cell_meta,
        var=pd.DataFrame(
            {"chrom": bin_chrom, "start": bin_start, "end": bin_end},
            index=pd.Index(bin_id, name="chrom5k"),
        ),
    )
    return adata


def adata_to_df(adata, var_dim, obs_dim="cell", dtype=np.uint8):
    # reduce dtype specifically for atac
    max_value = np.iinfo(dtype).max  # prevent overflow
    adata.X[adata.X > max_value] = max_value
    adata.X = adata.X.astype(dtype)

    df = adata.to_df()
    df.index.name = obs_dim
    df.columns.name = var_dim
    return df


def snap_to_xarray(snap_path, bin_sizes=(5000,), gene=True, dtype=np.uint8):
    total_records = {}
    for bin_size in bin_sizes:
        adata = read_snap(snap_path, bin_kind=bin_size)
        var_dim = f"chrom{bin_size}".replace("000", "k")
        data = adata_to_df(adata=adata, var_dim=var_dim, obs_dim="cell", dtype=dtype)
        total_records[f"{var_dim}_da"] = data
    if gene:
        adata = read_snap(snap_path, bin_kind="gene")
        data = adata_to_df(adata=adata, var_dim="gene", obs_dim="cell", dtype=dtype)
        total_records[f"gene_da"] = data
    ds = xr.Dataset(total_records)
    return ds


def snap_to_zarr(
    snap_path,
    output_path,
    bin_sizes=(5000,),
    gene=True,
    dtype=np.uint8,
    index_prefix=None,
):
    ds = snap_to_xarray(snap_path, bin_sizes=bin_sizes, gene=gene, dtype=dtype)
    if index_prefix is not None:
        # add a prefix to the index to make it unique when concatenate across different datasets
        ds.coords["cell"] = index_prefix + "_" + ds.get_index("cell")
    ds.to_zarr(output_path)
    return
