import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.multitest import multipletests


def _calculate_enrichment_score(raw_adata, labels):
    """
    Enrichment score modified from :cite:p:`Zeisel2018` for normalized methylation fractions
    Assuming the methylation value is posterior frac calculated by MCDS.add_mc_frac)
    """
    n_cells = raw_adata.shape[0]
    total_hypo_n = pd.Series(np.ravel((raw_adata.X < 1).sum(axis=0)).copy(),
                             index=raw_adata.var_names)
    total_sum = pd.Series(np.ravel(raw_adata.X.sum(axis=0).copy()),
                          index=raw_adata.var_names)

    enrichment = []
    label_orders = []
    for label, sub_df in raw_adata.obs.groupby(labels):
        in_cells = sub_df.shape[0]
        out_cells = n_cells - in_cells
        if in_cells < 2:
            continue

        label_orders.append(label)
        # hypo-methylated cells in the cluster
        in_hypo_n = pd.Series(np.ravel(
            (raw_adata[sub_df.index, :].X < 1).sum(axis=0)).copy(),
                              index=raw_adata.var_names)
        in_hypo_frac = in_hypo_n / in_cells
        # hypo-methylated cells not in the cluster
        out_hypo_n = total_hypo_n - in_hypo_n
        out_hypo_frac = out_hypo_n / out_cells

        # mean of cells in the cluster
        in_mean = pd.Series(np.ravel(
            raw_adata[sub_df.index, :].X.mean(axis=0)).copy(),
                            index=raw_adata.var_names)
        # mean of cells not in the cluster
        out_mean = (total_sum - in_mean * in_cells) / out_cells

        # finally, the enrichment score
        # In Zeisel et al., the frac and mean change on the same direction,
        # but in the mc frac case, larger fraction means smaller mean
        # so the fraction part is reversed
        e = ((in_hypo_frac + 0.1) / (out_hypo_frac + 0.1) *
             (out_mean + 0.01) / (in_mean + 0.01))
        enrichment.append(e)
        # enrichment > 1, the gene is hypo-methylated in that cluster, the larger the better

    enrichment = pd.DataFrame(enrichment,
                              index=label_orders).loc[:, raw_adata.var_names].copy()
    return enrichment


def _calculate_enrichment_score_cytograph(adata, labels):
    """
    The original CEF algorithm from :cite:p:`Zeisel2018` for count based data (RNA, ATAC)
    """
    n_cells = adata.shape[0]

    clusters = []
    nnz = []
    nnz_other = []
    means = []
    means_other = []
    sizes = []  # Number of cells per cluster
    total_nnz = np.ravel((adata.X > 0).sum(axis=0))
    total_sum = np.ravel(adata.X.sum(axis=0))
    total_cells = adata.shape[0]
    for cluster, sub_df in adata.obs.groupby(labels):
        clusters.append(cluster)
        this_cells = sub_df.shape[0]
        sizes.append(this_cells)
        is_in_cluster = adata.obs_names.isin(sub_df.index)
        # Number of nonzero values per cluster and other cells
        this_nnz = np.ravel((adata[is_in_cluster, :].X > 0).sum(axis=0))
        nnz.append(this_nnz)
        other_nnz = total_nnz - this_nnz
        nnz_other.append(other_nnz)
        # Mean value per cluster
        this_mean = np.ravel(adata[is_in_cluster, :].X.mean(axis=0))
        means.append(this_mean)
        other_mean = (total_sum - this_mean * this_cells) / \
                     (total_cells - this_cells)
        means_other.append(other_mean)
    nnz = np.vstack(nnz).astype(np.float64)
    nnz_other = np.vstack(nnz_other).astype(np.float64)
    means = np.vstack(means)
    means_other = np.vstack(means_other)
    sizes = np.array(sizes)

    # turn nnz count into fractions
    _sizes = sizes[:, None]
    _sizes_other = n_cells - _sizes
    nnz /= _sizes
    nnz_other /= _sizes_other

    enrichment = ((nnz + 0.1) / (nnz_other + 0.1) * (means + 0.01) /
                  (means_other + 0.01))
    enrichment = pd.DataFrame(enrichment,
                              index=clusters,
                              columns=adata.var_names)
    return enrichment


def _plot_enrichment_result(qvals, enrichment, null_enrichment, alpha):
    """Make some plots for the p-value, q-value and # of CEFs distribution"""
    fig, (f_ax, q_ax, e_ax) = plt.subplots(ncols=3, figsize=(15, 5))

    # plot feature with q_value pass cutoff
    ax = f_ax
    sns.histplot((qvals < alpha).sum(axis=0), ax=ax)
    ax.set_xlabel("# of clusters")
    ax.set_title(f"Features enriched in # of clusters (q<{alpha})")
    sns.despine(ax=ax)

    # plot q_value
    ax = q_ax
    sns.histplot(np.ravel(qvals), ax=ax, bins=50)
    ax.set(yticks=[],
           ylabel="",
           xlabel="q value",
           title="q value distribution")
    sns.despine(ax=ax, left=True)

    # plot enrichment distribution
    ax = e_ax
    plot_data = pd.DataFrame({
        "enrichment":
            enrichment.unstack().reset_index(drop=True),
        "null":
            null_enrichment.unstack().reset_index(drop=True),
    })
    plot_data = plot_data.unstack().reset_index()
    plot_data.columns = ["Type", "_", "Enrichment Score"]
    sns.histplot(
        data=plot_data,
        x="Enrichment Score",
        hue="Type",
        ax=ax,
        bins=100,
        binrange=enrichment.unstack().dropna().quantile((0.0, 0.95)).tolist(),
    )
    ax.set(yticks=[], ylabel="", title="Enrichment score distribution")
    sns.despine(ax=ax, left=True)

    fig.tight_layout()
    return


def _aggregate_enrichment(adata, enrichment, top_n, alpha, qvals, cluster_col):
    """Aggregate enrichment results, calculate q values"""
    # calculate use_features without qvals
    use_features_no_q = []
    for _, row in enrichment.iterrows():
        use_features_no_q.append(
            row.sort_values(ascending=False)[:top_n].dropna())
    use_features_no_q = pd.concat(use_features_no_q).index.unique()

    # mask non-sig enrichment scores, calculate again
    enrichment[qvals > alpha] = np.nan
    use_features = []
    for _, row in enrichment.iterrows():
        use_features.append(row.sort_values(ascending=False)[:top_n].dropna())
    use_features = pd.concat(use_features).index.unique()
    print(f"Selected {use_features.size} unique features")
    if len(use_features) == 0:
        print(f"No features found significantly enriched, "
              f"use top enriched features that did not pass q<{alpha}")
        use_features = use_features_no_q

    adata.uns[f"{cluster_col}_feature_enrichment"] = {
        "qvals": np.array(qvals).T,
        "enrichment": np.array(enrichment).T,
        "cluster_order": enrichment.index.tolist(),
    }
    return use_features


def cluster_enriched_features(adata: anndata.AnnData,
                              cluster_col: str,
                              top_n=200,
                              alpha=0.05,
                              stat_plot=True,
                              method="mc"):
    """
    Calculate top Cluster Enriched Features (CEF) from per-cell normalized dataset.
    An post-clustering feature selection step adapted from :cite:p:`Zeisel2018,La_Manno2021`
    and their great [cytograph2](https://github.com/linnarsson-lab/cytograph2) package.
    For details about CEF calculation, read the methods of :cite:p:`Zeisel2018`. Note that
    in original paper, they look for cluster-specific highly expressed genes as CEFs;
    for methylation, we are looking for hypo-methylation as CEFs, so the score and test is reversed.

    Parameters
    ----------
    adata
        adata containing per-cell normalized values.
        For methylation fraction, the value need to be 1-centered
        (1 means cell's average methylation), like those produced by
        :func:`ALLCools.mcds.mcds.MCDS.add_mc_frac` with `normalize_per_cell=True`.
        For RNA and ATAC, you can use per cell normalized counts.
        Do not log transform the data before running this function
    cluster_col
        The name of categorical variable in adata.obs
    top_n
        Select top N CEFs for each cluster
    alpha
        FDR corrected q-value cutoff
    stat_plot
        Whether making some summary plots for the CEF calculation
    method
        "mc" for methylation CEF (look for hypo-methylation),
        "rna" and "atac" for the RNA and ATAC or any count based data
        (use the original cytograph algorithm, look for higher value)

    Returns
    -------
    Modify adata inplace, adding a dictionary in adata.uns called f"{cluster_col}_feature_enrichment"
    The dictionary contains "qvals" (np.ndarray cluster-by-feature enrichment score q-value) and
    "cluster_order" (cluster order of the "qvals")
    """
    n_labels = adata.obs[cluster_col].unique().size
    print(f"Found {n_labels} clusters to compute feature enrichment score")

    if method == "mc":
        use_func = _calculate_enrichment_score
    elif method == "rna":
        use_func = _calculate_enrichment_score_cytograph
    elif method == "atac":
        use_func = _calculate_enrichment_score_cytograph
    else:
        use_func = _calculate_enrichment_score_cytograph
        print(
            f'method needs to be in ["mc", "rna", "atac"], got {method}. Using the algorithm for RNA.'
        )

    if n_labels == 1:
        print("No clusters detected, returning all features")
        use_features = adata.var_names.copy()
    else:
        print("Computing enrichment score")
        enrichment = use_func(adata, adata.obs[cluster_col])
        # shuffle cluster label to calculate null dist
        null_label = adata.obs[cluster_col].sample(n=adata.shape[0],
                                                   replace=False)
        null_label.index = adata.obs_names
        null_enrichment = use_func(adata, null_label)

        print("Computing enrichment score FDR-corrected P values")
        # FDR correction per row, row is cluster
        qvals = np.zeros_like(enrichment)
        for row in range(enrichment.shape[0]):
            null_values = null_enrichment.iloc[row, :].sort_values()
            values = enrichment.iloc[row, :]
            pvals = 1 - np.searchsorted(null_values, values) / values.shape[0]
            (_, q, _, _) = multipletests(pvals, alpha, method="fdr_bh")
            qvals[row, :] = q
        qvals = pd.DataFrame(qvals,
                             index=enrichment.index,
                             columns=enrichment.columns)

        if stat_plot:
            _plot_enrichment_result(
                qvals=qvals,
                enrichment=enrichment,
                null_enrichment=null_enrichment,
                alpha=alpha,
            )

        use_features = _aggregate_enrichment(
            adata=adata,
            enrichment=enrichment,
            top_n=top_n,
            alpha=alpha,
            qvals=qvals,
            cluster_col=cluster_col,
        )

    # save the calculated results in the input adata
    adata.var[f"{cluster_col}_enriched_features"] = adata.var_names.isin(
        use_features)
    return
