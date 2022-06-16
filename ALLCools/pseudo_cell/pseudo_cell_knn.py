import numpy as np
import warnings
from .pseudo_cell_kmeans import _merge_pseudo_cell


from enum import Enum
class SamplingStrategyEnum(Enum):
    SPARSE: str = 'sparse'
    DENSE: str = 'dense'

class ExamplerAndNeighborhoodSampler:
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
        '''
        Returns a general distance threshold for neighborhood of give size by random sampling.

                Parameters:
                        pulp_size (int): size of the neighborhood
                        n_tests (int): Sampling number
                        pulp_thicken_ratio (float): to get a more robust distance estimate, the neighborhood size 
                            will be relax by this ratio during sampling
                        robust_quantile (float): to get a more robust distance estimate, this quantile of the all
                            sampled distance will be used as the final distance
                Returns:
                        dist_thresh (str): the distance threshold
        '''        
        dists = self._sample_pulp_dist(n_tests, int(pulp_thicken_ratio * pulp_size))
        dist_thresh = np.quantile(dists, robust_quantile)
        return dist_thresh

            
    def _sample_fruit(
        self,
        pulp_size : int,
        n_kernels : int = None,
        max_attempts : int = 100,
        dist_thresh : float = None,
        ovlp_tol : float = 0.2,
        min_pulp_size : int = None,
        k=1000,
        strategy : SamplingStrategyEnum ='dense',
    ):
        '''
        Returns a list of examplers and their neighborhoods.
            Iteratively generating new examplers and corresponding neighborhoods. 
            In each iteration, at most one new exampler and the neighborhood is added to the final list.

                Parameters:
                        n_kernels (int): desired number of examplers; as many as posible if None
                        pulp_size (int): desired size of the neighborhood
                        max_attempts (int): max number of iterations of find a new exampler
                        dist_thresh (float): min distance allowed between examplers
                        ovlp_tol (float): max overlapping ratio between neighbors
                        min_pulp_size (int): for some examplers, the neighborhood sizes may not 
                            reach the desired size. This is the mininum acceptable neighobrhood size.
                        k (int): numbers of new exampler candidates to consider in each iteration. 
                            No need to change this argument in general.
                        strategy (['sparse','dense']): strategy of adding a new exampler.
                            'dense' will choose the candidate closest to the examplers already-selected. 
                                Recomand to combine 'dense' strategy with 'n_kernels = None'
                            'sparse' will choose the candidate fartherest to the examplers already-selected. 
                                Recomand to use 'sparse' strategy when a specific 'n_kernels' is desired.
                Returns:
                        kernels (list): indexes of selected examplers
                        pulps (list of list): list of indexes of the corresponding neighborhoods
        '''      
        
        if strategy==SamplingStrategyEnum.DENSE and n_kernels is not None:
            warnings.warn(
                    "'dense' strategy is better to use with 'n_kernels = None' "
                    "to get the best cover of the data"
                )
        
        import scipy.spatial.distance as ssd

        _dist_thresh = (
            dist_thresh if dist_thresh is not None else _select_dist_thresh(pulp_size)
        )

        kernels = []
        pulps = []
        unused = set(range(len(self.pca)))
        while (
            (max_attempts > 0) and 
            (n_kernels is None or len(kernels) < n_kernels)
        ):
            if len(kernels) > 0:
                kernel_cands = np.random.choice(list(unused), k)
                dists = ssd.cdist(self.pca[kernel_cands], self.pca[kernels])
                dists = dists.min(1)

                dists, kernel_cands = zip(*sorted(zip(dists, kernel_cands)))
                dists = np.array(dists)
                kernel_cands = np.array(kernel_cands)
                kernel_cands = kernel_cands[dists > _dist_thresh]
                dists = dists[dists > _dist_thresh]
                
                if len(dists)==0:
                    max_attempts-=1
                    continue
                    
                _, kernel_cands = zip(*sorted(zip(dists,kernel_cands), 
                                              reverse=(strategy==SamplingStrategyEnum.SPARSE)
                                             ))
            else:
                kernel_cands = [np.random.choice(list(unused))]

            changed = False
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
                        changed = True
                        break
            if not changed:
                max_attempts-=1

        return kernels, pulps

    def sample_examplers_and_neighborhoods(
        self,
        n_examplers : int,
        n_neighbors : int,
        min_n_neighbors : int = None,
        ovlp_tol : float = 0,
        dist_thresh : float = None,
        strategy : SamplingStrategyEnum = 'dense',
        max_attempts : int = 100,
    ):
        if dist_thresh is None:
            dist_thresh = self._select_dist_thresh(n_neighbors)
        if min_n_neighbors is None:
            min_n_neighbors = n_neighbors

        examplers, neighbors = self._sample_fruit(
            pulp_size=n_neighbors,
            n_kernels=n_examplers,
            max_attempts=max_attempts,
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
    n_pseudos=None,
    strategy : SamplingStrategyEnum = 'dense'
    
):
    _cell_meta = cell_meta[[cluster_col]].copy()
    index_name = _cell_meta.index.name
    _cell_meta = _cell_meta.reset_index()
    small_cluster_flags = []

    for c, cmeta in _cell_meta.groupby(cluster_col, as_index=False):
        _n_pseudos = len(cmeta) // target_pseudo_size
        if _n_pseudos == 0:
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
            sampler = ExamplerAndNeighborhoodSampler(coords[cmeta.index], n_components)
            pseudo_centers, pseudo_groups = sampler.sample_examplers_and_neighborhoods(
                n_pseudos, target_pseudo_size, min_pseudo_size, ovlp_tol=pseudo_ovlp, strategy='dense',
            )
        for i, (pcenter, pgroup) in enumerate(zip(pseudo_centers, pseudo_groups)):
            _cell_meta.loc[cmeta.iloc[pcenter].name, "pseudo_center"] = f"{c}::{i}"
#             _cell_meta.loc[cmeta.iloc[pgroup].index, "pseudo_cell"] = f"{c}::{i}"
            _cell_meta.loc[cmeta.iloc[pgroup].index, f"pseudo_cell::{c}::{i}"] = f"{c}::{i}"

    _cell_meta = _cell_meta.set_index(index_name)

    stats = _cell_meta.copy()
    stats["pseudo_cell"] = (~stats.loc[:, stats.columns.str.startswith("pseudo_cell::")].isna()).any()
    stats.index.name = "total_cells"
    stats = stats.reset_index().groupby(cluster_col, as_index=False).count()
    stats["cover_ratio"] = stats["pseudo_cell"] / stats["total_cells"]
    stats = stats.loc[:,~stats.columns.str.startswith("pseudo_cell::")]
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
    pseudo_cell_adata.obs[cluster_col] = pseudo_cell_adata.obs_names.str.split("::").str[0]
    return pseudo_cell_adata
