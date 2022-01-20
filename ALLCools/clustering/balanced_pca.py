"""
Perform Balanced PCA as well as
Reproducible PCA (A class allow user provide fitted model and just transform MCDS to adata with PCs)
"""
import anndata
import joblib
from sklearn.preprocessing import StandardScaler, RobustScaler
import scanpy as sc
import numpy as np
from scipy.stats import ks_2samp
import pandas as pd


def log_scale(
        adata, method="standard", with_mean=False, with_std=True, max_value=10, scaler=None
):
    """
    Perform log transform and then scale the cell-by-feature matrix

    Parameters
    ----------
    adata
        adata with normalized, unscaled cell-by-feature matrix
    method
        the type of scaler to use:
        'standard' for :class:`sklearn.preprocessing.StandardScaler`;
        'robust' for :class:`sklearn.preprocessing.RobustScaler`.
    with_mean
        Whether scale with mean center
    with_std
        Whether scale the std
    max_value
        Whether clip large values after scale
    scaler
        A fitted sklearn scaler, if provided, will only use it to transform the adata.

    Returns
    -------
    adata.X is scaled in place, the fitted scaler object will be return if the `scaler` parameter is None.
    """
    # log transform
    if "log" in adata.uns and adata.uns["log"]:
        # already log transformed
        print("adata.X is already log transformed, skip log step.")
    else:
        adata.X = np.log(adata.X)
        adata.uns["log"] = True

    if scaler is not None:
        print("Using user provided scaler.")
        if isinstance(scaler, str):
            import joblib

            scaler = joblib.load(scaler)
        user_provide = True
    else:
        if method == "robust":
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


def significant_pc_test(
        adata: anndata.AnnData,
        p_cutoff=0.1, update=True, obsm="X_pca", downsample=50000
):
    """
    Perform two-sample Kolmogorov-Smirnov test for goodness of fit on two adjacent PCs,
    select top PCs based on the `p_cutoff`. Top PCs have significantly different distributions, while
    later PCs only capturing random noise will have larger p-values. An idea from :cite:p:`Zeisel2018`.

    Parameters
    ----------
    adata
        adata with PC matrix calculated and stored in adata.obsm
    p_cutoff
        the p-value cutoff to select top PCs
    update
        Whether modify adata.obsm and only keep significant PCs
    obsm
        name of the PC matrix in adata.obsm
    downsample
        If the dataset is too large, downsample the cells before testing.

    Returns
    -------
    n_components
        number of PCs selected
    """
    pcs = adata.obsm[obsm]
    if pcs.shape[0] > downsample:
        print(
            f"Downsample PC matrix to {downsample} cells to calculate significant PC components"
        )
        use_pcs = pd.DataFrame(pcs).sample(downsample).values
    else:
        use_pcs = pcs
    i = 0
    for i in range(use_pcs.shape[1] - 1):
        cur_pc = use_pcs[:, i]
        next_pc = use_pcs[:, i + 1]
        _, p = ks_2samp(cur_pc, next_pc)
        if p > p_cutoff:
            break
    n_components = min(i + 1, use_pcs.shape[1])
    min_pc = min(4, use_pcs.shape[1])
    if n_components < min_pc:
        print(f'only {n_components} passed the P cutoff, '
              f'in order to proceed following analysis, will use first {min_pc} PCs')
        n_components = min_pc
    else:
        print(f"{n_components} components passed P cutoff of {p_cutoff}.")
    if update:
        adata.obsm[obsm] = pcs[:, :n_components]
        print(
            f"Changing adata.obsm['X_pca'] from shape {pcs.shape} to {adata.obsm[obsm].shape}"
        )
    return n_components


def balanced_pca(
        adata: anndata.AnnData, groups: str = "pre_clusters", max_cell_prop=0.1, n_comps=200, scale=False
):
    """
    Given a categorical variable (e.g., a pre-clustering label), perform balanced PCA by downsample
    cells in the large categories to make the overall population more balanced, so the PCs are expected
    to represent more variance among small categories.

    Parameters
    ----------
    adata
        adata after preprocessing and feature selection steps
    groups
        the name of the categorical variable in adata.obsm
    max_cell_prop
        any single category with cells > `n_cell * max_cell_prop` will be downsampled to this number.
    n_comps
        Number of components in PCA
    scale
        whether to scale the input matrix before PCA

    Returns
    -------
    adata with PC information stored in obsm, varm and uns like the :func:`scanpy.tl.pca` do.
    """
    # downsample large clusters
    use_cells = []
    size_to_downsample = max(int(adata.shape[0] * max_cell_prop), 50)
    for cluster, sub_df in adata.obs.groupby(groups):
        if sub_df.shape[0] > size_to_downsample:
            use_cells += sub_df.sample(
                size_to_downsample, random_state=0
            ).index.tolist()
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
    sc.tl.pca(
        adata_train,
        n_comps=n_comps,
        zero_center=True,
        svd_solver="arpack",
        random_state=0,
        return_info=False,
        use_highly_variable=None,
        dtype="float32",
        copy=False,
        chunked=False,
        chunk_size=None,
    )

    # transfer PCA result to full adata
    if downsample:
        if scale:
            adata.X = scaler.transform(adata.X)  # scale all cells with the same scaler
        adata.varm["PCs"] = adata_train.varm["PCs"]
        adata.obsm["X_pca"] = adata.X @ adata_train.varm["PCs"]
        adata.uns["pca"] = adata_train.uns["pca"]
    return adata


def get_pc_centers(adata: anndata.AnnData,
                   group: str,
                   outlier_label=None,
                   obsm="X_pca"):
    """
    Get the cluster centroids of the PC matrix

    Parameters
    ----------
    adata
        adata with cluster labels in obs and PC in obsm
    group
        the name of cluster labels in adata.obs
    outlier_label
        if there are outlier labels in obs, will exclude outliers before calculating centroids.
    obsm
        the key of PC matrix in obsm

    Returns
    -------
    pc_center
        a dataframe for cluster centroids by PC
    """
    pc_matrix = adata.obsm[obsm]
    pc_center = (
        pd.DataFrame(pc_matrix, index=adata.obs_names)
            .groupby(adata.obs[group])
            .median()
    )
    if outlier_label is not None:
        pc_center = pc_center[pc_center.index != outlier_label]
    return pc_center


class ReproduciblePCA:
    def __init__(
            self,
            scaler,
            mc_type,
            adata=None,
            pca_obj=None,
            pc_loading=None,
            var_names=None,
            max_value=10,
    ):
        """
        Perform preprocessing, feature selection and PCA using fitted scaler and PC transform

        Parameters
        ----------
        scaler :
            fitted sklearn scaler, can be the one return from
            {func}`log_scale <ALLCools.clustering.balanced_pca.log_scale>`
        mc_type :
            value to select mc_type dim in MCDS
        adata :
            If anndata is provided, will take the PC loading from adata, otherwise, pca_obj or pc_loading is needed
        pca_obj :
            fitted sklearn PCA model, will take the pc_loading from the pca_obj.components_
        pc_loading :
            or you can directly provide the pc_loading matrix
        var_names :
            the var_names of the pc_loading matrix, will use these index to select features
        max_value :
            the maximum to clip after scaling
        """
        if adata is not None:
            self.pc_loading = adata.varm["PCs"]
            self.pc_vars = adata.var_names
        else:
            if (pca_obj is None) and (pc_loading is None):
                raise ValueError("Need to provide either PCA object or PC Loadings")
            if pca_obj is not None:
                self.pc_loading = pca_obj.components_
            else:
                self.pc_loading = pc_loading
            if var_names is None:
                raise ValueError("Need to provide var_names")
            self.pc_vars = var_names

        self.scaler = scaler
        self.var_dim = self.pc_vars.name
        self.mc_type = mc_type
        self.max_value = max_value

    def mcds_to_adata(self, mcds):
        """
        Get adata from MCDS with only selected features

        Parameters
        ----------
        mcds : Input raw count MCDS object

        Returns
        -------
        adata:
            Adata with per-cell normalized mC fraction and selected features
        """
        if f"{self.var_dim}_da_frac" not in mcds:
            # adding mC frac and normalize per cell should be done separately for every cell/dataset
            mcds.add_mc_frac(
                var_dim=self.var_dim, normalize_per_cell=True, clip_norm_value=10
            )

        # add hvf info into mcds
        hvf_judge = mcds.get_index(self.var_dim).isin(self.pc_vars)
        mcds.coords[f"{self.var_dim}_{self.mc_type}_feature_select"] = hvf_judge

        # adata only contains self.pc_vars
        adata = mcds.get_adata(
            mc_type=self.mc_type,
            var_dim=self.var_dim,
            da_suffix="frac",
            obs_dim="cell",
            select_hvf=True,
            split_large_chunks=True,
        )
        return adata

    def scale(self, adata):
        """
        Perform {func}`log_scale <ALLCools.clustering.balanced_pca.log_scale>` with fitted scaler

        Parameters
        ----------
        adata :
            Adata with per-cell normalized mC fraction and selected features

        Returns
        -------
        adata.X is transformed in place
        """
        log_scale(
            adata,
            method="standard",
            with_mean=False,
            with_std=True,
            max_value=self.max_value,
            scaler=self.scaler,
        )

    def pc_transform(self, adata):
        """
        calculate the PC from adata.X and PC loading, store PCs in adata.obsm["X_pca"] and
        loadings in adata.varm["PCs"]

        Parameters
        ----------
        adata :
            Adata with log_scale transformed mC fraction and selected features
        Returns
        -------
        PC information stored in adata.obsm and varm
        """
        pc = np.dot(adata.X, self.pc_loading)
        adata.obsm["X_pca"] = pc
        adata.varm["PCs"] = self.pc_loading

    def mcds_to_adata_with_pc(self, mcds):
        """
        From raw count MCDS to adata object with PCs using fitted scaler and PC loadings.
        Steps include select features, calculate per-cell normalized mC fractions,
        log_scale transform the data with fitted scaler, and finally add PC matrix.

        Parameters
        ----------
        mcds :
            Raw count MCDS

        Returns
        -------
        adata :
            anndata object with per-cell normalized, log scale transformed matrix in .X and PCs in
            adata.obsm["X_pca"] and PC loadings in adata.varm["PCs"]. The scale and PC are done with fitted model.
        """
        adata = self.mcds_to_adata(mcds)
        self.scale(adata)
        self.pc_transform(adata)
        return adata

    def dump(self, path):
        """
        Save the ReproduciblePCA to path
        """
        joblib.dump(self, path)
