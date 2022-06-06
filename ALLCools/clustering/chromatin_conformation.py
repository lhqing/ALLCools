import seaborn as sns


def filter_feature_by_dispersion(adata, top_n=0.5, plot=True, inplace=False):
    """
    Filter features based on dispersion.
    """
    n_features = int(adata.shape[1] * top_n)
    n_features = max(n_features, 1)

    mean = adata.X.mean(axis=0)
    mean[mean == 0] += 1e-10
    dispersion = adata.X.var(axis=0) / mean
    adata.var['dispersion'] = dispersion
    use_features = adata.var['dispersion'].sort_values(ascending=False)[:n_features].index
    adata.var['dispersion_filtered'] = adata.var_names.isin(use_features)

    if plot:
        low, high = adata.var['dispersion'].quantile([0.01, 0.99])
        delta = (high - low) / 10
        bin_range = (low - delta, high + delta)
        sns.displot(data=adata.var,
                    x='dispersion',
                    bins=100,
                    binrange=bin_range,
                    hue='dispersion_filtered')

    if inplace:
        adata._inplace_subset_var(adata.var['dispersion_filtered'])
    return
