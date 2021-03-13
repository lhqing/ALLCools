import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl


def plot_decomp_scatters(adata,
                         n_components,
                         base_name='PC',
                         obsm_name='X_pca',
                         hue=None,
                         palette='viridis',
                         hue_quantile=(0.25, 0.75),
                         nrows=5,
                         ncols=5):
    available_comps = adata.obsm[obsm_name].shape[1]
    nrows = min(nrows, available_comps // 2 // ncols + 1)
    fig, axes = plt.subplots(figsize=(ncols * 3, nrows * 3),
                             nrows=nrows,
                             ncols=ncols,
                             dpi=150)
    vmin = None
    vmax = None
    for i, ax in enumerate(axes.ravel()):
        _x = i * 2
        _y = i * 2 + 1
        if _y > available_comps - 1:
            ax.axis('off')
            continue
        _plot_data = pd.DataFrame(
            {
                f'{base_name}{_x + 1}': adata.obsm[obsm_name][:, _x],
                f'{base_name}{_y + 1}': adata.obsm[obsm_name][:, _y]
            },
            index=adata.obs_names)
        if hue is not None:
            _plot_data[hue] = adata.obs[hue]
            hue_dtype = str(_plot_data[hue].dtype)
            if hue_dtype.startswith('float') or hue_dtype.startswith('int'):
                vmin, vmax = _plot_data[hue].quantile(hue_quantile).tolist()
                hue_norm = (vmin, vmax)
            else:
                hue_norm = None
        else:
            hue_norm = None
        sns.scatterplot(ax=ax,
                        data=_plot_data,
                        x=f'{base_name}{_x + 1}',
                        y=f'{base_name}{_y + 1}',
                        hue=hue,
                        hue_norm=hue_norm,
                        s=1,
                        linewidth=0,
                        palette=palette,
                        legend=None)

        # adjust axis
        xmin, xmax = _plot_data[f'{base_name}{_x + 1}'].quantile([0.01, 0.99
                                                                  ]).tolist()
        delta = (xmax - xmin) * 0.2
        xmin -= delta
        xmax += delta
        ymin, ymax = _plot_data[f'{base_name}{_y + 1}'].quantile([0.01, 0.99
                                                                  ]).tolist()
        delta = (ymax - ymin) * 0.2
        ymin -= delta
        ymax += delta
        ax.set(xticks=[], yticks=[], xlim=(xmin, xmax), ylim=(ymin, ymax))
        if _x + 1 <= n_components:
            ax.set_xlabel(f'{base_name}{_x + 1}', color='red')
        if _y + 1 <= n_components:
            ax.set_ylabel(f'{base_name}{_y + 1}', color='red')
    print(f'Red axis labels are used {base_name}s')
    if vmin is not None and vmax is not None:
        cnorm = mpl.colors.Normalize(vmin, vmax)
        fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm, cmap='viridis'),
                     ax=fig.axes,
                     shrink=0.6,
                     fraction=0.06,
                     label=hue)
    return fig, axes