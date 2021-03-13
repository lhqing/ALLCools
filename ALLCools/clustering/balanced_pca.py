from sklearn.preprocessing import StandardScaler
import scanpy as sc
import numpy as np
from scipy.stats import ks_2samp
import pandas as pd


def log_scale(adata):
    adata.X = np.log(adata.X)
    sc.pp.scale(adata, zero_center=False)
    return


def significant_pc_test(adata, p_cutoff=0.1, update=True, obsm='X_pca'):
    """

    Parameters
    ----------
    adata
    p_cutoff
    update

    Returns
    -------

    """
    pcs = adata.obsm[obsm]

    i = 0
    for i in range(pcs.shape[1] - 1):
        cur_pc = pcs[:, i]
        next_pc = pcs[:, i + 1]
        p = ks_2samp(cur_pc, next_pc).pvalue
        if p > p_cutoff:
            break
    n_components = min(i + 1, pcs.shape[1])
    print(f'{n_components} components passed P cutoff of {p_cutoff}.')
    if update:
        adata.obsm[obsm] = pcs[:, :n_components]
        print(f"Changing adata.obsm['X_pca'] from shape {pcs.shape} to {adata.obsm[obsm].shape}")
    return n_components


def balanced_pca(adata, groups='pre_clusters', max_cell_prop=0.1, n_comps=200):
    # downsample large clusters
    use_cells = []
    size_to_downsample = max(int(adata.shape[0] * max_cell_prop), 30)
    for cluster, sub_df in adata.obs.groupby(groups):
        if sub_df.shape[0] > size_to_downsample:
            use_cells += sub_df.sample(size_to_downsample, random_state=0).index.tolist()
        else:
            use_cells += sub_df.index.tolist()

    # get training adata
    if len(use_cells) == adata.shape[0]:
        downsample = False
        adata_train = adata
    else:
        downsample = True
        adata_train = adata[use_cells, :].copy()

    # pca
    scaler = StandardScaler()
    adata_train.X = scaler.fit_transform(adata_train.X)
    sc.pp.scale(adata_train)  # only scale training cells
    sc.tl.pca(adata_train,
              n_comps=n_comps,
              zero_center=True,
              svd_solver='arpack',
              random_state=0,
              return_info=False,
              use_highly_variable=None,
              dtype='float32',
              copy=False,
              chunked=False,
              chunk_size=None)

    # transfer PCA result to full adata
    if downsample:
        adata.X = scaler.transform(adata.X)  # scale all cells with the same scaler
        adata.varm['PCs'] = adata_train.varm['PCs']
        adata.obsm['X_pca'] = adata.X @ adata_train.varm['PCs']
        adata.uns['pca'] = adata_train.uns['pca']
    return adata


def get_pc_centers(adata, group, outlier_label=None, obsm='X_pca'):
    pc_matrix = adata.obsm[obsm]
    pc_center = pd.DataFrame(pc_matrix, index=adata.obs_names).groupby(
        adata.obs[group]).median()
    if outlier_label is not None:
        pc_center = pc_center[pc_center.index != outlier_label]
    return pc_center
