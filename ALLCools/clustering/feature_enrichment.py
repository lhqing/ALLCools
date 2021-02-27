import scanpy as sc
from scipy.stats import ks_2samp
import numpy as np
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def _simple_leiden(adata, k, n_components, resolution, plot):
    sc.pp.scale(adata)

    sc.tl.pca(adata,
              n_comps=100,
              zero_center=True,
              svd_solver='arpack',
              random_state=0,
              return_info=False,
              use_highly_variable=None,
              dtype='float32',
              copy=False,
              chunked=False,
              chunk_size=None)

    if n_components == 'auto':
        # determine automatically
        for i in range(adata.obsm['X_pca'].shape[1] - 1):
            cur_pc = adata.obsm['X_pca'][:, i]
            next_pc = adata.obsm['X_pca'][:, i + 1]
            p = ks_2samp(cur_pc, next_pc).pvalue
            if p > 0.1:
                break
        n_components = i + 1
    else:
        n_components = int(n_components)
        assert n_components > 1
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=n_components)
    sc.tl.leiden(adata, resolution=resolution)

    if plot:
        sc.tl.umap(adata)
        sc.pl.umap(adata, color='leiden')
    return adata


def calculate_enrichment_score(raw_adata, labels):
    """
    Enrichment score modified from Zeisel et al. 2018 cell
    """
    n_cells = raw_adata.shape[0]
    sizes = labels.value_counts()

    # hypo-methylated cells in the cluster
    in_hypo_n = pd.DataFrame(raw_adata.X < 1,
                             index=raw_adata.obs_names,
                             columns=raw_adata.var_names).groupby(labels).sum()
    in_hypo_frac = in_hypo_n / sizes.values[:, None]
    # hypo-methylated cells not in the cluster
    freature_sum = in_hypo_n.sum(axis=0)
    out_hypo_n = freature_sum.values[None, :] - in_hypo_n
    out_sizes = n_cells - sizes
    out_hypo_frac = out_hypo_n / out_sizes.values[:, None]

    # mean of cells in the cluster
    in_mean = pd.DataFrame(raw_adata.X,
                           index=raw_adata.obs_names,
                           columns=raw_adata.var_names).groupby(labels).mean()
    # mean of cells not in the cluster
    out_mean = (raw_adata.X.sum(axis=0)[None, :] -
                in_mean * sizes.values[:, None]) / out_sizes.values[:, None]

    # finally, the enrichment score
    enrichment = (in_hypo_frac + 0.1) / (out_hypo_frac + 0.1) * (
            in_mean + 0.01) / (out_mean + 0.01)
    return enrichment


def _plot_enrichment_result(qvals, enrichment, null_enrichment, alpha):
    fig, (f_ax, q_ax, e_ax) = plt.subplots(ncols=3, figsize=(15, 5))

    # plot feature with q_value pass cutoff
    ax = f_ax
    sns.histplot((qvals < alpha).sum(axis=0), ax=ax)
    ax.set_xlabel('# of clusters')
    ax.set_title(f'Features enriched in # of clusters (q<{alpha})')
    sns.despine(ax=ax)

    # plot q_value
    ax = q_ax
    sns.histplot(qvals.ravel(), ax=ax)
    ax.set(yticks=[], ylabel='', xlabel='q value', title='q value distribution')
    sns.despine(ax=ax, left=True)

    # plot enrichment distribution
    ax = e_ax
    plot_data = pd.DataFrame({
        'enrichment':
            enrichment.unstack().reset_index(drop=True),
        'null':
            null_enrichment.unstack().reset_index(drop=True)
    })
    plot_data = plot_data.unstack().reset_index()
    plot_data.columns = ['Type', '_', 'Enrichment Score']
    sns.histplot(data=plot_data, x='Enrichment Score', hue='Type', ax=ax)
    ax.set(yticks=[], ylabel='', title='Enrichment score distribution')
    sns.despine(ax=ax, left=True)

    fig.tight_layout()
    return


def feature_enrichment(adata, top_n=100, k=25, n_components='auto', resolution=1, alpha=0.05, plot=True,
                       cluster_plot=False):
    raw_adata = adata
    # computation happens on a copy of adata, input adata is unchanged
    print('Computing clusters')
    adata = _simple_leiden(adata.copy(), k=k, n_components=n_components, resolution=resolution, plot=cluster_plot)
    n_labels = adata.obs['leiden'].unique().size
    print(f'Found {n_labels} clusters to compute feature enrichment score')
    if n_labels == 1:
        print('No clusters detected, returning all features')
        use_features = adata.var_names.copy()

    else:
        print('Computing enrichment score')
        enrichment = calculate_enrichment_score(raw_adata, adata.obs['leiden'])
        null_label = adata.obs.sample(n=adata.shape[0])['leiden']
        null_label.index = adata.obs_names
        null_enrichment = calculate_enrichment_score(raw_adata, null_label)

        print('Computing enrichment score FDR-corrected P values')
        qvals = np.zeros_like(enrichment)
        for ix in range(enrichment.shape[1]):
            null_values = null_enrichment.iloc[:, ix].sort_values()
            values = enrichment.iloc[:, ix]
            pvals = 1 - np.searchsorted(null_values, values) / values.shape[0]
            (_, q, _, _) = multipletests(pvals, alpha, method="fdr_bh")
            qvals[:, ix] = q

        if plot:
            _plot_enrichment_result(qvals=qvals, enrichment=enrichment, null_enrichment=null_enrichment, alpha=alpha)

        # calculate use_features without qvals
        use_features_no_q = []
        for _, row in enrichment.iterrows():
            use_features_no_q.append(row.sort_values(ascending=False)[:top_n].dropna())
        use_features_no_q = pd.concat(use_features_no_q).index.unique()

        # mask non-sig enrichment scores, calculate again
        enrichment[qvals > 0.05] = np.nan
        use_features = []
        for _, row in enrichment.iterrows():
            use_features.append(row.sort_values(ascending=False)[:top_n].dropna())
        use_features = pd.concat(use_features).index.unique()
        print(f'Selected {use_features.size} unique features')
        if len(use_features) == 0:
            print(f'No features found significantly enriched, '
                  f'use top enriched features that did not pass q<{alpha}')
            use_features = use_features_no_q

    raw_adata.obs['pre_clusters'] = adata.obs['leiden']
    raw_adata.var['pre_clusters_enriched'] = adata.var_names.isin(use_features)
    return
