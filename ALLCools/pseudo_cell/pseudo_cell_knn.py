import numpy as np
import warnings
from .pseudo_cell_kmeans import _merge_pseudo_cell


class ContractedExamplerSampler:
    def __init__(self, data, n_components=30, normalize=False):
        self.data = data
        if self.data.shape[1] > n_components:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler

            data_std = StandardScaler().fit_transform(self.data)
            self.pca = PCA(n_components).fit_transform(data_std)
        else:
            self.pca = np.array(data)

        if normalize:
            # from sklearn.preprocessing import MaxAbsScaler
            # self.pca = MaxAbsScaler().fit_transform(self.pca)
            raise NotImplementedError

        from pynndescent import NNDescent

        self.ann_index = NNDescent(self.pca)

    def _sample_pulp_dist(self, n_kernels, pulp_size):
        kernels = np.random.choice(len(self.pca), n_kernels)
        pulps, dists = self.ann_index.query(self.pca[kernels], pulp_size)
        return dists

    def _select_dist_thresh(
        self, pulp_size, n_tests=100, pulp_thicken_ratio=1.2, robust_quantile=0.9
    ):
        dists = self._sample_pulp_dist(n_tests, int(pulp_thicken_ratio * pulp_size))
        return np.quantile(dists, robust_quantile)

    def _sample_fruit(
        self,
        n_kernels,
        pulp_size,
        max_iters,
        dist_thresh=None,
        ovlp_tol=0.2,
        min_pulp_size=None,
        k=100,
    ):
        import scipy.spatial.distance as ssd

        _dist_thresh = (
            dist_thresh if dist_thresh is not None else _select_dist_thresh(pulp_size)
        )

        kernels = []
        pulps = []
        unused = set(range(len(self.pca)))
        while max_iters > 0 and len(kernels) < n_kernels:
            max_iters -= 1
            if len(kernels) > 0:
                kernel_cands = np.random.choice(list(unused), k)
                dists = ssd.cdist(self.pca[kernel_cands], self.pca[kernels])
                dists = dists.min(1)

                dists, kernel_cands = zip(*sorted(zip(dists, kernel_cands)))
                kernel_cands = np.array(kernel_cands)
                # dists = np.array(dists)
                kernel_cands = kernel_cands[dists > _dist_thresh]
            else:
                kernel_cands = [np.random.choice(list(unused))]

            for kernel in kernel_cands:
                pulp, dists = self.ann_index.query([self.pca[kernel]], pulp_size)
                pulp = pulp.flatten()

                if (
                    (dist_thresh is None)
                    or ((min_pulp_size is None) and (dists < dist_thresh).all())
                    or (
                        (min_pulp_size is not None)
                        and (dists < dist_thresh).sum() >= min_pulp_size
                    )
                ):

                    if len(set(pulp) - set(unused)) / len(pulp) <= ovlp_tol:
                        kernels.append(kernel)
                        pulps.append(pulp)
                        unused -= set(pulp)
                        break

        return kernels, pulps

    def sample_contracted_examplers(
        self,
        n_examplers,
        n_neighbors,
        min_n_neighbors=None,
        ovlp_tol=0,
        dist_thresh=None,
        max_iters=100,
    ):
        if dist_thresh is None:
            dist_thresh = self._select_dist_thresh(n_neighbors)

        examplers, neighbors = self._sample_fruit(
            n_kernels=n_examplers,
            pulp_size=n_neighbors,
            max_iters=max_iters,
            dist_thresh=dist_thresh,
            ovlp_tol=ovlp_tol,
            min_pulp_size=min_n_neighbors,
            k=100,
        )
        return examplers, neighbors


def sample_pseudo_cells(
    cell_meta,
    cluster_col,
    coords,
    target_pseudo_size,
    min_pseudo_size=None,
    ignore_small_cluster=False,
    n_components=30,
    pseudo_ovlp=0,
):
    _cell_meta = cell_meta[[cluster_col]].copy()
    index_name = _cell_meta.index.name
    _cell_meta = _cell_meta.reset_index()
    small_cluster_flags = []

    for c, cmeta in _cell_meta.groupby(cluster_col, as_index=False):
        n_pseudos = len(cmeta) // target_pseudo_size
        if n_pseudos == 0:
            if ignore_small_cluster:
                continue
            else:
                warnings.warn(
                    f'Size of cluster "{c}" is smaller than target pseudo-cell size.'
                )
                small_cluster_flags.append(True)
                pseudo_centers = [0]
                pseudo_groups = [list(range(cmeta.shape[0]))]
        else:
            small_cluster_flags.append(False)
            sampler = ContractedExamplerSampler(coords[cmeta.index], n_components)
            pseudo_centers, pseudo_groups = sampler.sample_contracted_examplers(
                n_pseudos, target_pseudo_size, min_pseudo_size, ovlp_tol=pseudo_ovlp
            )
        for i, (pcenter, pgroup) in enumerate(zip(pseudo_centers, pseudo_groups)):
            _cell_meta.loc[cmeta.iloc[pcenter].name, "pseudo_center"] = f"{c}::{i}"
            _cell_meta.loc[cmeta.iloc[pgroup].index, "pseudo_cell"] = f"{c}::{i}"

    _cell_meta = _cell_meta.set_index(index_name)

    stats = _cell_meta.copy()
    stats.index.name = "total_cells"
    stats = stats.reset_index().groupby(cluster_col, as_index=False).count()
    stats["cover_ratio"] = stats["pseudo_cell"] / stats["total_cells"]
    stats.columns = [
        cluster_col,
        "total_cells",
        "pseudo_cells",
        "covered_cells",
        "cover_ratio",
    ]

    # stats.index = ['total_cells', 'pseudo_cells', 'covered_cells']
    # stats['pseud_yield'] = stats['covered_cells']/stats['total_cells']

    return _cell_meta, stats


# a wrapper of sample_pseudo_cells to use adata as input
def generate_pseudo_cells(
    adata,
    cluster_col="leiden",
    obsm="X_pca",
    target_pseudo_size=100,
    min_pseudo_size=None,
    ignore_small_cluster=False,
    n_components=None,
    aggregate_func="downsample",
    pseudo_ovlp=0,
):
    if n_components is None:
        n_components = adata.obsm[obsm].shape[1]
    if min_pseudo_size is None:
        min_pseudo_size = 1

    # determine cell group
    cell_group, stats = sample_pseudo_cells(
        cell_meta=adata.obs,
        cluster_col=cluster_col,
        coords=adata.obsm[obsm],
        target_pseudo_size=target_pseudo_size,
        min_pseudo_size=min_pseudo_size,
        ignore_small_cluster=ignore_small_cluster,
        n_components=n_components,
        pseudo_ovlp=pseudo_ovlp,
    )
    adata.obs["cell_group"] = cell_group["pseudo_cell"]
    pseudo_cell_adata = _merge_pseudo_cell(adata=adata, aggregate_func=aggregate_func)
    pseudo_cell_adata.obs[cluster_col] = pseudo_cell_adata.obs_names.str.split(
        "::"
    ).str[0]
    return pseudo_cell_adata
