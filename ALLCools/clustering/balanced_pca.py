import joblib
from sklearn.preprocessing import StandardScaler, RobustScaler
import scanpy as sc
import numpy as np
from scipy.stats import ks_2samp
import pandas as pd


def log_scale(adata, method='standard', with_mean=False, with_std=True, max_value=10, scaler=None):
    # log transform
    if 'log' in adata.uns and adata.uns['log']:
        # already log transformed
        print('adata.X is already log transformed, skip log step.')
    else:
        adata.X = np.log(adata.X)
        adata.uns['log'] = True

    if scaler is not None:
        print('Using user provided scaler.')
        if isinstance(scaler, str):
            import joblib
            scaler = joblib.load(scaler)
        user_provide = True
    else:
        if method == 'robust':
            scaler = RobustScaler(with_centering=with_mean, with_scaling=with_std)
        else:
            scaler = StandardScaler(with_mean=with_mean, with_std=with_std)
        user_provide = False

    # transform data
    if user_provide:
        adata.X = scaler.transform(adata.X)
    else:
        adata.X = scaler.fit_transform(adata.X)

    # clip large values
    if max_value is not None:
        adata.X[adata.X > max_value] = max_value
        adata.X[adata.X < -max_value] = -max_value

    # return scaler or not
    if user_provide:
        return
    else:
        return scaler


def significant_pc_test(adata, p_cutoff=0.1, update=True, obsm='X_pca', downsample=50000):
    """

    Parameters
    ----------
    adata
    p_cutoff
    update
    obsm
    downsample

    Returns
    -------

    """
    pcs = adata.obsm[obsm]
    if pcs.shape[0] > downsample:
        print(f'Downsample PC matrix to {downsample} cells to calculate significant PC components')
        use_pcs = pd.DataFrame(pcs).sample(downsample).values
    else:
        use_pcs = pcs
    i = 0
    for i in range(use_pcs.shape[1] - 1):
        cur_pc = use_pcs[:, i]
        next_pc = use_pcs[:, i + 1]
        p = ks_2samp(cur_pc, next_pc).pvalue
        if p > p_cutoff:
            break
    n_components = min(i + 1, use_pcs.shape[1])
    print(f'{n_components} components passed P cutoff of {p_cutoff}.')
    if update:
        adata.obsm[obsm] = pcs[:, :n_components]
        print(f"Changing adata.obsm['X_pca'] from shape {pcs.shape} to {adata.obsm[obsm].shape}")
    return n_components


def balanced_pca(adata, groups='pre_clusters', max_cell_prop=0.1, n_comps=200, scale=False):
    # downsample large clusters
    use_cells = []
    size_to_downsample = max(int(adata.shape[0] * max_cell_prop), 50)
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

    # in case cells are smaller than n_comps
    n_comps = min(min(adata_train.shape), n_comps)

    # scale (optional)
    if scale:
        scaler = StandardScaler()
        adata_train.X = scaler.fit_transform(adata_train.X)
    else:
        scaler = None

    # pca
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
        if scale:
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


class ReproduciblePCA:
    def __init__(self, scaler, mc_type, adata=None, pca_obj=None, pc_loading=None, var_names=None, max_value=10):
        if adata is not None:
            self.pc_loading = adata.varm['PCs']
            self.pc_vars = adata.var_names
        else:
            if (pca_obj is None) and (pc_loading is None):
                raise ValueError('Need to provide either PCA object or PC Loadings')
            if pca_obj is not None:
                self.pc_loading = pca_obj.components_
            else:
                self.pc_loading = pc_loading
            if var_names is None:
                raise ValueError('Need to provide var_names')
            self.pc_vars = var_names

        self.scaler = scaler
        self.var_dim = self.pc_vars.name
        self.mc_type = mc_type
        self.max_value = max_value

    def mcds_to_adata(self, mcds):
        if f'{self.var_dim}_da_frac' not in mcds:
            # adding mC frac and normalize per cell should be done separately for every cell/dataset
            mcds.add_mc_frac(var_dim=self.var_dim,
                             normalize_per_cell=True,
                             clip_norm_value=10)

        # add hvf info into mcds
        hvf_judge = mcds.get_index(self.var_dim).isin(self.pc_vars)
        mcds.coords[f'{self.var_dim}_{self.mc_type}_feature_select'] = hvf_judge

        adata = mcds.get_adata(mc_type=self.mc_type,
                               var_dim=self.var_dim,
                               da_suffix='frac',
                               obs_dim='cell',
                               select_hvf=True,
                               split_large_chunks=True)
        return adata

    def scale(self, adata):
        log_scale(adata,
                  method='standard',
                  with_mean=False,
                  with_std=True,
                  max_value=self.max_value,
                  scaler=self.scaler)

    def pc_transform(self, adata):
        pc = np.dot(adata.X, self.pc_loading)
        adata.obsm['X_pca'] = pc
        adata.varm['PCs'] = self.pc_loading

    def mcds_to_adata_with_pc(self, mcds):
        adata = self.mcds_to_adata(mcds)
        self.scale(adata)
        self.pc_transform(adata)
        return adata

    def dump(self, path):
        joblib.dump(self, path)
