import numpy as np
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def calculate_enrichment_score(raw_adata, labels):
    """
    Enrichment score modified from Zeisel et al. 2018 Cell for normalized methylation fractions
    """
    n_cells = raw_adata.shape[0]
    sizes = labels.value_counts()

    # hypo-methylated cells in the cluster
    in_hypo_n = pd.DataFrame(raw_adata.X < 1,
                             index=raw_adata.obs_names,
                             columns=raw_adata.var_names).groupby(labels).sum()
    in_hypo_frac = in_hypo_n / sizes.values[:, None]
    # hypo-methylated cells not in the cluster
    feature_sum = in_hypo_n.sum(axis=0)
    out_hypo_n = feature_sum.values[None, :] - in_hypo_n
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
    # In Zeisel et al, the frac and mean change on the same direction,
    # but in the mc frac case, larger fraction means smaller mean
    # so the fraction part is reversed
    enrichment = ((out_hypo_frac + 0.1) / (in_hypo_frac + 0.1)) * \
                 ((in_mean + 0.01) / (out_mean + 0.01))
    # enrichment direction is the same as mC fraction
    # enrichment < 1, the gene is hypo-methylated in that cluster
    # enrichment > 1, the gene is hyper-methylated in that cluster
    return enrichment


def calculate_enrichment_score_cytograph(adata, labels):
    """
    Enrichment score algorithm from Zeisel et al. 2018 Cell
    """
    n_cells = adata.shape[0]

    # Number of cells per cluster
    sizes = labels.value_counts()

    clusters = []
    nnz = []
    means = []
    for cluster, sub_df in adata.obs.groupby(labels):
        clusters.append(cluster)
        # Number of nonzero values per cluster
        nnz.append(np.ravel((adata[sub_df.index].X > 0).sum(axis=0)))
        # Mean value per cluster
        means.append(np.ravel(adata[sub_df.index].X.mean(axis=0)))
    nnz = np.vstack(nnz)
    means = np.vstack(means)

    nnz_overall = nnz.sum(axis=0)
    means_overall = np.ravel(adata.X.mean(axis=0))

    _sizes = sizes.values[:, None]
    f_nnz = nnz / _sizes
    f_nnz_overall = nnz_overall / n_cells
    # Means and fraction non-zero values in other clusters (per cluster)
    means_other = ((means_overall * n_cells)[None] -
                   (means * _sizes)) / (n_cells - _sizes)
    means_other[means_other < 0] = 0
    f_nnz_other = ((f_nnz_overall * n_cells)[None] -
                   (f_nnz * _sizes)) / (n_cells - _sizes)
    f_nnz_other[f_nnz_other < 0] = 0

    enrichment = (f_nnz + 0.1) / (f_nnz_other +
                                  0.1) * (means + 0.01) / (means_other + 0.01)
    enrichment = pd.DataFrame(enrichment, index=clusters, columns=adata.var_names)
    return enrichment


def _plot_enrichment_result(qvals, enrichment, null_enrichment, alpha):
    fig, (f_ax, q_ax, e_ax) = plt.subplots(ncols=3, figsize=(15, 5))

    # plot feature with q_value pass cutoff
    ax = f_ax
    sns.histplot((qvals < alpha).sum(axis=0), bins=30, ax=ax)
    ax.set_xlabel('# of clusters')
    ax.set_title(f'Features enriched in # of clusters (q<{alpha})')
    sns.despine(ax=ax)

    # plot q_value
    ax = q_ax
    sns.histplot(qvals.ravel(), ax=ax, bins=50)
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
    sns.histplot(data=plot_data,
                 x='Enrichment Score',
                 hue='Type',
                 ax=ax,
                 bins=100,
                 binrange=enrichment.unstack().dropna().quantile(
                     (0., 0.95)).tolist())
    ax.set(yticks=[], ylabel='', title='Enrichment score distribution')
    sns.despine(ax=ax, left=True)

    fig.tight_layout()
    return


def _aggregate_enrichment(adata, enrichment, top_n, alpha, qvals, cluster_col):
    # calculate use_features without qvals
    use_features_no_q = []
    for _, row in enrichment.iterrows():
        use_features_no_q.append(row.sort_values(ascending=False)[:top_n].dropna())
    use_features_no_q = pd.concat(use_features_no_q).index.unique()

    # mask non-sig enrichment scores, calculate again
    enrichment[qvals > alpha] = np.nan
    use_features = []
    for _, row in enrichment.iterrows():
        use_features.append(row.sort_values(ascending=False)[:top_n].dropna())
    use_features = pd.concat(use_features).index.unique()
    print(f'Selected {use_features.size} unique features')
    if len(use_features) == 0:
        print(f'No features found significantly enriched, '
              f'use top enriched features that did not pass q<{alpha}')
        use_features = use_features_no_q

    adata.uns[f'{cluster_col}_feature_enrichment'] = {'qvals': qvals.T,
                                                      'cluster_order': enrichment.index.tolist()}
    return use_features


def cluster_enriched_features(adata,
                              cluster_col,
                              top_n=200,
                              alpha=0.05,
                              stat_plot=True,
                              method='mc'):
    """

    Parameters
    ----------
    adata
    cluster_col
    top_n
    alpha
    stat_plot
    method
        "mc" for methylation fraction,
        "rna" and "atac" for the original algorithm for the RNA

    Returns
    -------

    """
    n_labels = adata.obs[cluster_col].unique().size
    print(f'Found {n_labels} clusters to compute feature enrichment score')

    if method == 'mc':
        use_func = calculate_enrichment_score
    elif method == 'rna':
        use_func = calculate_enrichment_score_cytograph
    elif method == 'atac':
        use_func = calculate_enrichment_score_cytograph
    else:
        use_func = calculate_enrichment_score_cytograph
        print(f'method needs to be in ["mc", "rna", "atac"], got {method}. Using the algorithm for RNA.')

    if n_labels == 1:
        print('No clusters detected, returning all features')
        use_features = adata.var_names.copy()
    else:
        print('Computing enrichment score')
        enrichment = use_func(adata, adata.obs[cluster_col])
        null_label = adata.obs.sample(n=adata.shape[0])[cluster_col]
        null_label.index = adata.obs_names
        null_enrichment = use_func(adata, null_label)

        print('Computing enrichment score FDR-corrected P values')
        qvals = np.zeros_like(enrichment)
        for ix in range(enrichment.shape[1]):
            null_values = null_enrichment.iloc[:, ix].sort_values()
            values = enrichment.iloc[:, ix]
            pvals = 1 - np.searchsorted(null_values, values) / values.shape[0]
            (_, q, _, _) = multipletests(pvals, alpha, method="fdr_bh")
            qvals[:, ix] = q

        if stat_plot:
            _plot_enrichment_result(qvals=qvals, enrichment=enrichment, null_enrichment=null_enrichment, alpha=alpha)

        use_features = _aggregate_enrichment(adata=adata,
                                             enrichment=enrichment,
                                             top_n=top_n,
                                             alpha=alpha,
                                             qvals=qvals,
                                             cluster_col=cluster_col)
    # save the calculated results in the input adata
    adata.var[f'{cluster_col}_enriched_features'] = adata.var_names.isin(use_features)
    return
