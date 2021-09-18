import numpy as np
import warnings


class ContractedExamplerSampler:
    def __init__(self, data, n_components=30, normalize=False):
        self.data = data
        if self.data.shape[1]>n_components:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            data_std = StandardScaler().fit_transform(self.data)
            self.pca = PCA(n_components).fit_transform(data_std)
        else:
            self.pca = np.array(data)
            
        if normalize:
#             from sklearn.preprocessing import MaxAbsScaler
#             self.pca = MaxAbsScaler().fit_transform(self.pca)
            raise NotImplementedError
        
        from pynndescent import NNDescent
        self.ann_index = NNDescent(self.pca)

        
    def _sample_pulp_dist(self, n_kernels, pulp_size):
        kernels = np.random.choice(len(self.pca), n_kernels)
        pulps, dists = self.ann_index.query(self.pca[kernels], pulp_size)
        return dists
        
        
    def _select_dist_thresh(self, pulp_size, n_tests=100, pulp_thicken_ratio=1.2, robust_quantile=0.9):
        dists = self._sample_pulp_dist(n_tests, int(pulp_thicken_ratio*pulp_size))
        return np.quantile(dists, robust_quantile)
    
    
    def _sample_fruit(self, n_kernels, pulp_size, max_iters, dist_thresh=None, 
                      ovlp_tol=0.2, min_pulp_size=None, k=100):
        
        import scipy.spatial.distance as ssd
        
        _dist_thresh = dist_thresh if dist_thresh is not None else _select_dist_thresh(pulp_size)
        
        kernels = []
        pulps = []
        
        unused = set(range(len(self.pca)))
        
        while max_iters>0 and len(kernels)<n_kernels:
            max_iters-=1

            if len(kernels)>0:
                kernel_cands = np.random.choice(list(unused),k)
                dists = ssd.cdist(self.pca[kernel_cands], self.pca[kernels])
                dists = dists.min(1)
                
                dists, kernel_cands = zip(*sorted(zip(dists, kernel_cands)))
                kernel_cands = np.array(kernel_cands)
#                 dists = np.array(dists)
                kernel_cands = kernel_cands[dists>_dist_thresh]
            else:
                kernel_cands = [np.random.choice(list(unused))]
                
                
            for kernel in kernel_cands:
                pulp, dists = self.ann_index.query([self.pca[kernel]], pulp_size)
                pulp = pulp.flatten()


                if (dist_thresh is None) or \
                    ((min_pulp_size is None) and (dists<dist_thresh).all()) or \
                    (dists<dist_thresh).sum()>=min_pulp_size:

                    if len(set(pulp)-set(unused))/len(pulp)<= ovlp_tol:
                        kernels.append(kernel)
                        pulps.append(pulp)
                        unused-=set(pulp)
                        break
                        
        return kernels, pulps
    
    
    def sample_contracted_examplers(self, n_examplers, n_neighbors, min_n_neighbors=None, ovlp_tol=0, 
                                    dist_thresh=None, max_iters=100, ):
        if dist_thresh is None:
            dist_thresh = self._select_dist_thresh(n_neighbors)
        
        examplers, neighbors = self._sample_fruit(n_kernels=n_examplers, pulp_size=n_neighbors, 
                                                  max_iters=max_iters, dist_thresh=dist_thresh, 
                                                  ovlp_tol=ovlp_tol, min_pulp_size=min_n_neighbors, 
                                                  k=100)
        return examplers, neighbors


def sample_pseudo_cells(cell_meta, cluster_col, coords, tgt_pseudo_size, min_pseudo_size=None, 
                        ignore_small_cluster=False, n_components=30):

    _cell_meta = cell_meta[[cluster_col]].copy()
    index_name = _cell_meta.index.name
    _cell_meta = _cell_meta.reset_index()
    small_cluster_flags = []
    
    for c, cmeta in _cell_meta.groupby(cluster_col, as_index=False):
        n_pseudos = len(cmeta)//tgt_pseudo_size
        if n_pseudos == 0:
            if ignore_small_cluster:
                continue
            else:
                warnings.warn(f'Size of cluster "{c}" is smaller than target psuedo cell size.')
                small_cluster_flags.append(True)
                pseudo_centers = [cmeta.index.tolist()[0]]
                pseudo_groups = [cmeta.index.tolist()]
        else:
            small_cluster_flags.append(False)
            sampler = ContractedExamplerSampler(coords[cmeta.index], n_components)
            pseudo_centers, pseudo_groups = sampler.sample_contracted_examplers(n_pseudos, tgt_pseudo_size, 
                                                                                min_pseudo_size, ovlp_tol=0)
        for i, (pcenter, pgroup) in enumerate(zip(pseudo_centers, pseudo_groups)):
            _cell_meta.loc[cmeta.loc[pcenter].name, 'pseudo_center'] = f'{c}::{i}'
            _cell_meta.loc[cmeta.loc[pgroup].index, 'pseudo_cell'] = f'{c}::{i}'
        
    _cell_meta = _cell_meta.set_index(index_name)

    
    stats = _cell_meta.copy()
    stats.index.name = 'total_cells'
    stats = stats.reset_index().groupby(cluster_col, as_index=False).count()
    stats['cover_ratio'] = stats['pseudo_cell']/stats['total_cells']
    stats.columns = [cluster_col, 'total_cells','pseudo_cells','covered_cells','cover_ratio']

    #stats.index = ['total_cells', 'pseudo_cells', 'covered_cells']
    #stats['pseud_yield'] = stats['covered_cells']/stats['total_cells']
    
    return _cell_meta, stats
    
