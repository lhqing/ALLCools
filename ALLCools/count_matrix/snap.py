import h5py
import anndata
import pandas as pd
import numpy as np
import scipy.sparse as ss

SNAP_BD_FIELDS = {
    'TN': 'TotalFrag',
    'UM': 'UniqueMapFrag',
    'SE': 'SingleEndFrag',
    'SA': 'SecondaryAlign',
    'PE': 'PairEndFrag',
    'PP': 'ProperlyPairedFrag',
    'PL': 'ProperFragLength',
    'US': 'UsableFrag',
    'UQ': 'UniqFrag',
    'CM': 'chrMFrag'
}


def read_snap(file_path, bin_kind=5000):
    """Read snap hdf5 file into anndata.AnnData. """
    if bin_kind == 'gene':
        return _read_snap_genes(file_path)
    elif isinstance(bin_kind, int):
        return _read_snap_bins(file_path, bin_size=bin_kind)
    else:
        raise NotImplementedError


def _read_snap_meta(f):
    cell_barcode = np.array(f['/BD/name']).astype(str)
    cell_meta = pd.DataFrame([], index=cell_barcode)
    for name in f['/BD/'].keys():
        if name == 'name':
            continue
        if name in SNAP_BD_FIELDS:
            col_name = SNAP_BD_FIELDS[name]
        else:
            col_name = name
        cell_meta[col_name] = np.array(f[f'/BD/{name}'])
    return cell_meta, cell_barcode


def _read_snap_genes(file_path):
    """Read gene data from snap hdf5 file into anndata.AnnData"""

    with h5py.File(file_path) as f:
        data = np.array(f[f'/GM/count'])
        idx = np.array(f[f'/GM/idx']) - 1  # R is 1 based, python is 0 based
        idy = np.array(f[f'/GM/idy']) - 1  # R is 1 based, python is 0 based
        gene_id = np.array(f[f'/GM/name']).astype(str)

        cell_meta, cell_barcode = _read_snap_meta(f)
        data = ss.coo_matrix((data, (idx, idy)), shape=(cell_barcode.size, gene_id.size)).tocsc()

    adata = anndata.AnnData(X=data,
                            obs=cell_meta,
                            var=pd.DataFrame(index=pd.Index(gene_id, name='gene')))
    return adata


def _read_snap_bins(file_path, bin_size=5000):
    """Read bin data from snap hdf5 file into anndata.AnnData"""
    with h5py.File(file_path) as f:
        data = np.array(f[f'/AM/{bin_size}/count'])
        idx = np.array(f[f'/AM/{bin_size}/idx']) - 1  # R is 1 based, python is 0 based
        idy = np.array(f[f'/AM/{bin_size}/idy']) - 1  # R is 1 based, python is 0 based

        bin_chrom = np.array(f[f'/AM/{bin_size}/binChrom']).astype(str)
        bin_start = np.array(f[f'/AM/{bin_size}/binStart']) - 1  # 0 based bed format
        bin_end = np.array(f[f'/AM/{bin_size}/binStart']) - 1 + bin_size  # 0 based bed format
        bin_id = np.core.defchararray.add(np.core.defchararray.add(bin_chrom, '-'),
                                          bin_start.astype(str))

        cell_meta, cell_barcode = _read_snap_meta(f)
        data = ss.coo_matrix((data, (idx, idy)), shape=(cell_barcode.size, bin_id.size)).tocsc()

    adata = anndata.AnnData(X=data,
                            obs=cell_meta,
                            var=pd.DataFrame({'chrom': bin_chrom,
                                              'start': bin_start,
                                              'end': bin_end},
                                             index=pd.Index(bin_id, name='chrom5k')))
    return adata
