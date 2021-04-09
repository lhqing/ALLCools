import pathlib
import anndata
import pandas as pd
import xarray as xr
import warnings
import scanpy as sc
import numpy as np
from sklearn.metrics import roc_auc_score
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from ..plot import categorical_scatter, continuous_scatter


def single_pairwise_dmg(cluster_l, cluster_r,
                        top_n, adj_p_cutoff,
                        delta_rate_cutoff,
                        auroc_cutoff,
                        adata_dir,
                        dmg_dir):
    # load data
    adata_l = anndata.read_h5ad(
        f'{adata_dir}/{cluster_l}.h5ad')
    adata_r = anndata.read_h5ad(
        f'{adata_dir}/{cluster_r}.h5ad')

    # generate single adata for DMG
    adata = adata_l.concatenate(adata_r,
                                batch_key='groups',
                                index_unique=None)
    adata.obs = pd.concat([adata_l.obs, adata_r.obs])
    try:
        assert adata.obs_names.duplicated().sum() == 0
    except AssertionError as e:
        print(cluster_l, cluster_r)
        raise e
    # reverse_adata, centered by 1 because after normalization all prior is center to 1
    adata.X = (adata.X - 1) * -1 + 1

    # calc DMG
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata,
                                groupby='groups',
                                n_genes=top_n,
                                method='wilcoxon')
    dmg_result = pd.DataFrame({
        data_key: pd.DataFrame(adata.uns['rank_genes_groups'][data_key],
                               columns=[cluster_r, cluster_l]).stack()
        for data_key in ['names', 'pvals_adj']
    })
    dmg_result.reset_index(drop=True, inplace=True)

    # annotate cluster_delta
    dmg_result['left-right'] = f'{cluster_l}-{cluster_r}'
    l_mean = pd.Series(np.mean(adata_l.X, axis=0), index=adata_l.var_names)
    left_mean = dmg_result['names'].map(l_mean)
    r_mean = pd.Series(np.mean(adata_r.X, axis=0), index=adata_l.var_names)
    right_mean = dmg_result['names'].map(r_mean)
    dmg_result['delta'] = left_mean - right_mean

    # filter
    dmg_result = dmg_result[(dmg_result['pvals_adj'] < adj_p_cutoff) & (
            dmg_result['delta'].abs() > delta_rate_cutoff)].copy()
    dmg_result['hypo_in'] = dmg_result['delta'].apply(
        lambda i: cluster_l if i < 0 else cluster_r)
    dmg_result['hyper_in'] = dmg_result['delta'].apply(
        lambda i: cluster_r if i < 0 else cluster_l)
    dmg_result = dmg_result.set_index('names').drop_duplicates()

    # add AUROC and filter again
    auroc = {}
    for gene, row in dmg_result.iterrows():
        yscore = adata.obs_vector(gene)
        ylabel = adata.obs['groups'] == row['hypo_in']
        score = roc_auc_score(ylabel, yscore)
        score = abs(score - 0.5) + 0.5
        auroc[gene] = score
    dmg_result['AUROC'] = pd.Series(auroc)
    dmg_result = dmg_result[(dmg_result['AUROC'] > auroc_cutoff)].copy()

    # save
    dmg_result.to_hdf(f'{dmg_dir}/{cluster_l}-{cluster_r}.hdf', key='data')
    return


class PairwiseDMG:
    def __init__(self,
                 max_cell_per_group=1000,
                 top_n=10000,
                 adj_p_cutoff=0.001,
                 delta_rate_cutoff=0.3,
                 auroc_cutoff=0.9,
                 random_state=0,
                 n_jobs=1):
        self.X = None
        self.groups = None
        self._obs_dim = ''
        self._var_dim = ''
        self.dmg_table = None
        self._outlier_label = None

        # parameters
        self.max_cell_per_group = max_cell_per_group
        self.top_n = top_n
        self.adj_p_cutoff = adj_p_cutoff
        self.delta_rate_cutoff = delta_rate_cutoff
        self.auroc_cutoff = auroc_cutoff
        self.random_state = random_state
        self.n_jobs = n_jobs

        # internal
        self._adata_dir = '_adata_for_dmg'
        self._dmg_dir = '_dmg_results'

    def fit_predict(self, x, groups, obs_dim='cell', var_dim='gene', outlier='Outlier', cleanup=True):
        if (len(x.shape) != 2) or not isinstance(x, xr.DataArray):
            raise ValueError('Expect an cell-by-feature 2D xr.DataArray as input matrix.')
        self._obs_dim = obs_dim
        self._var_dim = var_dim
        self._outlier_label = outlier

        self.X = x
        use_order = x.get_index(obs_dim)
        self.groups = groups.astype('category').reindex(use_order)

        # save adata for each group to dict
        print('Generating cluster AnnData files')
        self._save_cluster_adata()

        # run pairwise DMG
        print('Computing pairwise DMG')
        self._pairwise_dmg()

        # cleanup
        if cleanup:
            self._cleanup()

    def _save_cluster_adata(self):
        adata_dir = pathlib.Path(self._adata_dir)
        adata_dir.mkdir(exist_ok=True)

        for cluster, sub_series in self.groups.groupby(self.groups):
            if cluster == self._outlier_label:
                # skip outlier
                continue
            output_path = adata_dir / f'{cluster}.h5ad'
            if output_path.exists():
                continue
            if sub_series.size > self.max_cell_per_group:
                sub_series = sub_series.sample(self.max_cell_per_group,
                                               random_state=self.random_state)
            cluster_adata = anndata.AnnData(X=self.X.sel({self._obs_dim: sub_series.index}).values,
                                            obs=pd.DataFrame({'groups': sub_series.astype('category')}),
                                            var=pd.DataFrame([], index=self.X.get_index(self._var_dim)))
            cluster_adata.write_h5ad(output_path)
        return

    def _pairwise_dmg(self):
        dmg_dir = pathlib.Path(self._dmg_dir)
        dmg_dir.mkdir(exist_ok=True)

        pairs = [
            i for i in combinations(sorted(self.groups.unique()), 2)
            if self._outlier_label not in i
        ]
        print(len(pairs), 'pairwise DMGs')
        n_pairs = len(pairs)

        with ProcessPoolExecutor(min(n_pairs, self.n_jobs)) as exe:
            step = max(1, n_pairs // 10)
            futures = {}
            for (cluster_l, cluster_r) in pairs:
                f = exe.submit(single_pairwise_dmg,
                               cluster_l=cluster_l,
                               cluster_r=cluster_r,
                               top_n=self.top_n,
                               adj_p_cutoff=self.adj_p_cutoff,
                               delta_rate_cutoff=self.delta_rate_cutoff,
                               auroc_cutoff=self.auroc_cutoff,
                               adata_dir=self._adata_dir,
                               dmg_dir=self._dmg_dir)
                futures[f] = (cluster_l, cluster_r)
            for i, f in enumerate(as_completed(futures)):
                f.result()
                if i % step == 0:
                    print(f'{i + 1}/{n_pairs} finished')

        # summarize
        self.dmg_table = pd.concat((pd.read_hdf(p) for p in dmg_dir.glob('*.hdf')))

    def _cleanup(self):
        subprocess.run(['rm', '-rf', str(self._adata_dir), str(self._dmg_dir)], check=True)

    def aggregate_pairwise_dmg(self, adata, groupby, obsm='X_pca'):
        # using the cluster centroids in PC space to calculate dendrogram
        pc_matrix = adata.obsm[obsm]
        pc_center = pd.DataFrame(pc_matrix, index=adata.obs_names).groupby(adata.obs[groupby]).median()
        # calculate cluster pairwise similarity
        cluster_dist = pairwise_distances(pc_center)
        cluster_dist = pd.DataFrame(cluster_dist, index=pc_center.index, columns=pc_center.index)
        cluster_dist_norm = cluster_dist / cluster_dist.values.max()
        cluster_sim = 1 - cluster_dist_norm
        cluster_pair_sim_dict = {f'{a}-{b}': value for (a, b), value in cluster_sim.unstack().iteritems()}
        self.dmg_table['similarity'] = self.dmg_table['left-right'].map(cluster_pair_sim_dict)

        # aggregate pairwise DMG to get the cluster level DMG, use the similarity to normalize AUROC
        cluster_dmgs = {}
        for cluster, sub_df in self.dmg_table.groupby('hypo_in'):
            values = sub_df['AUROC'] * sub_df['similarity']
            cluster_dmg_order = values.groupby(values.index).sum().sort_values(ascending=False)
            cluster_dmgs[cluster] = cluster_dmg_order
        return cluster_dmgs


def single_ovr_dmg(cell_label, mcds, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff, auroc_cutoff):
    # get adata
    cell_judge = mcds.get_index(obs_dim).isin(cell_label.index)
    adata = mcds.sel({obs_dim: cell_judge}).get_adata(mc_type=mc_type, var_dim=var_dim, select_hvf=False)
    adata.var = pd.DataFrame([], index=adata.var_names)
    adata.obs['groups'] = cell_label.astype('category')

    # calc DMG
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata,
                                groupby='groups',
                                n_genes=top_n,
                                method='wilcoxon')
    dmg_result = pd.DataFrame({
        data_key: pd.DataFrame(adata.uns['rank_genes_groups'][data_key]).stack()
        for data_key in ['names', 'pvals_adj']
    })
    dmg_result = dmg_result[dmg_result.index.get_level_values(1).astype(bool)].reset_index(drop=True)
    # add fold change
    in_cells_mean = adata.X[adata.obs['groups'].astype(bool),].mean(axis=0)
    out_cells_mean = adata.X[~adata.obs['groups'].astype(bool),].mean(axis=0)
    fc = pd.Series(in_cells_mean / out_cells_mean, index=adata.var_names)
    dmg_result['fc'] = dmg_result['names'].map(fc)
    # filter
    dmg_result = dmg_result[(dmg_result['pvals_adj'] < adj_p_cutoff) & (
            dmg_result['fc'] < fc_cutoff)].copy()
    dmg_result = dmg_result.set_index('names').drop_duplicates()

    # add AUROC and filter again
    auroc = {}
    for gene, row in dmg_result.iterrows():
        yscore = adata.obs_vector(gene)
        ylabel = adata.obs['groups'] == True
        score = roc_auc_score(ylabel, yscore)
        score = abs(score - 0.5) + 0.5
        auroc[gene] = score
    dmg_result['AUROC'] = pd.Series(auroc)
    dmg_result = dmg_result[(dmg_result['AUROC'] > auroc_cutoff)].copy()
    return dmg_result


def one_vs_rest_dmg(cell_meta, group, mcds,
                    obs_dim='cell', var_dim='gene', mc_type='CHN',
                    top_n=1000, adj_p_cutoff=0.01, fc_cutoff=0.8, auroc_cutoff=0.8,
                    max_cluster_cells=2000, max_other_fold=5):
    """

    Parameters
    ----------
    cell_meta
    group
    mcds
    obs_dim
    var_dim
    mc_type
    top_n
    adj_p_cutoff
    fc_cutoff
    auroc_cutoff
    max_cluster_cells
    max_other_fold

    Returns
    -------

    """
    clusters = cell_meta[group].unique()
    dmg_table = []
    for cluster in sorted(clusters):
        print(f'Calculating cluster {cluster} DMGs.')
        # determine cells to use
        cluster_judge = cell_meta[group] == cluster
        in_cells = cluster_judge[cluster_judge]
        out_cells = cluster_judge[~cluster_judge]
        if in_cells.size > max_cluster_cells:
            in_cells = in_cells.sample(max_cluster_cells, random_state=0)

        max_other_cells = in_cells.size * max_other_fold
        if out_cells.size > max_other_cells:
            out_cells = out_cells.sample(max_other_cells, random_state=0)

        cell_label = pd.concat([in_cells, out_cells])
        dmg_df = single_ovr_dmg(cell_label, mcds, obs_dim, var_dim, mc_type, top_n, adj_p_cutoff, fc_cutoff,
                                auroc_cutoff)
        dmg_df['cluster'] = cluster
        dmg_table.append(dmg_df)
    dmg_table = pd.concat(dmg_table)
    return dmg_table
