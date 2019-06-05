import anndata
import scipy.sparse as ss
from ..utilities import parse_chrom_size, generate_chrom_bin_bed_dataframe


def _transform_single_h5ad(adata_path, output_path, chrom_size_path,
                           bin_size, step_size, window_size, compression):
    if (step_size % bin_size != 0) or (window_size % bin_size != 0):
        raise ValueError('step_size and window_size need to be integral multiple of bin_size')
    n = step_size // bin_size
    m = window_size // bin_size

    adata = anndata.read_h5ad(adata_path)

    # somehow, I need to copy this out otherwise its super slow
    chrom_idx = adata.var['chrom'].values.copy()
    csc_data = adata.X.tocsc()
    chrom_dict = parse_chrom_size(chrom_size_path)

    chrom_data_list = []
    for chrom in chrom_dict.keys():
        print(chrom)
        chrom_csc_data = csc_data[:, chrom_idx == chrom]
        chunk_generator = (ss.csc_matrix(chrom_csc_data[:, n:n + m].sum(axis=1))
                           for i in range(0, chrom_csc_data.shape[1], n))
        chrom_data = ss.hstack(list(chunk_generator))
        chrom_data_list.append(chrom_data)
    total_data = ss.hstack(chrom_data_list)

    # TODO add all necessary info in adata.uns
    adata = anndata.AnnData(X=total_data,
                            obs=adata.obs,
                            var=generate_chrom_bin_bed_dataframe(chrom_size_path,
                                                                 window_size=window_size,
                                                                 step_size=step_size),
                            uns=dict(bin_size=window_size,
                                     step_size=step_size,
                                     chrom_size_path=chrom_size_path))
    
    adata.write(filename=output_path, compression=compression)
    return output_path


# TODO add parallel version to transform all adata chunks and mc and covs
