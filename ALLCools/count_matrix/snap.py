import h5py
import anndata
import pandas as pd
import numpy as np
import scipy.sparse as ss


def read_snap(file_path, bin_size=5000):
    """Read snap hdf5 file into anndata.AnnData"""
    with h5py.File(file_path) as f:
        data = np.array(f[f'/AM/{bin_size}/count'])
        idx = np.array(f[f'/AM/{bin_size}/idx']) - 1  # R is 1 based, python is 0 based
        idy = np.array(f[f'/AM/{bin_size}/idy']) - 1  # R is 1 based, python is 0 based

        bin_chrom = np.array(f[f'/AM/{bin_size}/binChrom']).astype(str)
        bin_start = np.array(f[f'/AM/{bin_size}/binStart']) - 1  # 0 based bed format
        bin_end = np.array(f[f'/AM/{bin_size}/binStart']) - 1 + bin_size  # 0 based bed format
        bin_id = np.core.defchararray.add(np.core.defchararray.add(bin_chrom, '-'),
                                          bin_start.astype(str))

        cell_barcode = np.array(f['/BD/name']).astype(str)
        cell_id = [f'ATAC_{i}' for i in range(cell_barcode.size)]
        cell_meta = pd.DataFrame([], index=cell_id)
        cell_meta['barcode'] = cell_barcode
        for name in f['/BD/'].keys():
            if name == 'name':
                continue
            cell_meta[name] = np.array(f[f'/BD/{name}'])
        data = ss.coo_matrix((data, (idx, idy)), shape=(cell_barcode.size, bin_id.size)).tocsc()

    adata = anndata.AnnData(X=data,
                            obs=cell_meta,
                            var=pd.DataFrame({'chrom': bin_chrom,
                                              'start': bin_start,
                                              'end': bin_end},
                                             index=pd.Index(bin_id, name='chrom5k')))
    return adata
