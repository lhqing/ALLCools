import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from .dmg import PairwiseDMG
class ClusterMerge:
    def __init__(self, merge_criterion, stop_criterion=None, 
                 stop_clusters=-1, n_cells=200, metric='euclidean',method='average', 
                 label_concat_str='::'):

        self.data_for_tree = None
        self.cell_to_type = None
        self.gene_mcds = None

        
        self.n_cells = n_cells
        self.metric = metric
        self.method = method
        self.label_concat_str = label_concat_str
                
        self.merge_criterion = merge_criterion
        self.stop_criterion = stop_criterion
        self.stop_clusters = stop_clusters
        
        self.curr_tree = None
        self.curr_mean = None        
        self.merge_evidences = {}
        
    def _construct_tree(self):
        cells = self.cell_to_type.to_frame().groupby(self.cell_to_type.name)\
                    .apply(lambda x: x.sample(self.n_cells, replace=True)).droplevel(0).sort_index().index
        self.curr_mean = self.data_for_tree.loc[cells].groupby(self.cell_to_type[cells]).mean()
        pdist = scipy.spatial.distance.pdist(self.curr_mean, metric=self.metric)
        Z = sch.linkage(pdist, metric=self.metric, method=self.method)
        self.curr_tree = sch.to_tree(Z)
        
    @staticmethod
    def _traverse(node, call_back):
        if node.is_leaf():
            return
        elif node.left.is_leaf() and node.right.is_leaf():
            call_back(node)
        if not node.left.is_leaf():
            ClusterMerge._traverse(node.left, call_back)
        if not node.right.is_leaf():
            ClusterMerge._traverse(node.right, call_back)

    def _merge_pair(self, pair, concat_str='::'):
        left_id, right_id = pair
        left_lbl, right_lbl = self.curr_mean.iloc[[left_id, right_id]].index
        print('checking',left_lbl, right_lbl)#TODO
        pair_cells = self.cell_to_type[self.cell_to_type.isin([left_lbl, right_lbl])].index
        pair_cell_types = self.cell_to_type.loc[pair_cells]

        pair_mcds = self.gene_mcds.sel(cell = pair_cells)
        pair_mcds.load()

        separable, evidence, *_ = self.merge_criterion.predict((left_lbl, right_lbl), 
                                                               pair_cells, pair_cell_types, pair_mcds)
        mergeable = not separable


        if mergeable:
            self.cell_to_type.loc[pair_cells] = f'{left_lbl}{concat_str}{right_lbl}'
            
#             self.merge_evidences[left_lbl, right_lbl] = evidence
#             self.merge_evidences[right_lbl, left_lbl] = evidence
            self.merge_evidences[tuple(sorted([left_lbl, right_lbl]))] = evidence
    
            print(left_lbl, right_lbl, 'merged')#TODO

        return mergeable            
    
    
    def fit_predict(self, data_for_tree, cell_to_type, gene_mcds):
        
        self.data_for_tree = data_for_tree
        self.cell_to_type = cell_to_type
        self.gene_mcds = gene_mcds
        
        count = 0
        while True:
            count += 1
            self._construct_tree()
            pairs = []
            ClusterMerge._traverse(self.curr_tree, lambda x: pairs.append((x.left.id,x.right.id)))

            rlt = list(map(self._merge_pair, pairs))
            print('round',count, 'merged', sum(rlt))
            
            if self.stop_clusters>0 and len(self.cell_to_type.unique())<=self.stop_clusters:
                print()#TODO
                break
            if sum(rlt)==0:
                print()#TODO
                break
            if self.stop_criterion is not None \
                    and self.stop_criterion(self.data_for_tree, self.cell_to_type, self.gene_mcds):
                print()#TODO
                break
                
        return self.cell_to_type, self.merge_evidences


class PairwiseDMGCriterion:
    def __init__(self, 
                 max_cell_per_group = 100,
                 top_n_markers = 5,
                 adj_p_cutoff = 0.001,
                 delta_rate_cutoff = 0.3,
                 auroc_cutoff = 0.85,
                 use_modality = 'either',
                 random_state = 0,
                 n_jobs = 10, 
                 verbose = False,
                ):
        
        self.agg = {'either':np.logical_or,'both':np.logical_and}[use_modality]
        
        self.pwdmg = PairwiseDMG(max_cell_per_group = max_cell_per_group,
                                 top_n=top_n_markers,
                                 adj_p_cutoff=adj_p_cutoff,
                                 delta_rate_cutoff=delta_rate_cutoff,
                                 auroc_cutoff=auroc_cutoff,
                                 random_state=0,
                                 n_jobs=n_jobs,
                                 verbose=False,)        
#         self.max_cell_per_group = max_cell_per_group
        self.top_n_markers = top_n_markers
#         self.adj_p_cutoff = adj_p_cutoff
#         self.delta_rate_cutoff = delta_rate_cutoff
#         self.auroc_cutoff = auroc_cutoff
        self.use_modality = use_modality
#         self.n_jobs = n_jobs
        
        
    def predict(self, pair_labels, pair_cells, pair_cell_types, pair_mcds, da_name = 'gene_da_frac'):
        mc_types = ['CHN','CGN']
        separable = {x:False for x in mc_types}
        evidence = {}
        for mc_type in mc_types:
            self.pwdmg.fit_predict(x=pair_mcds[da_name].sel(mc_type=mc_type), 
                                   groups=pair_cell_types)
            evidence[mc_type] = self.pwdmg.dmg_table
            if len(self.pwdmg.dmg_table)>=self.top_n_markers:
                separable[mc_type] = True
                if self.use_modality=='either':
                    break
        return self.agg(*separable.values()), evidence
            

