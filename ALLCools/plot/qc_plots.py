import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import ticker


def plot_on_plate(data, value_col, groupby, ncols=4,
                  plate_base=384, figsize=(5, 3.5),
                  row_base='Row384',
                  col_base='Col384',
                  vmin=0, vmax=1, heatmap_kws=None, aggregation_func=None):
    """
    Plot metadata into 384 or 96 plate view (heatmap)
    Parameters
    ----------
    data
        dataframe contain all kinds of metadata
    value_col
        value to be plotted on plate view
    groupby
        groupby column, typically groupby plate id column(s) to plot each plate separately
    ncols
        number of column for axes, nrows will be calculated accordingly
    plate_base
        {384, 96} size of the plate view
    figsize
        matplotlib.Figure figsize
    vmin
        cmap vmin
    vmax
        cmap vmax
    heatmap_kws
        kws pass to sns.heatmap
    aggregation_func
        apply to reduce rows after groupby if the row is not unique
    """

    if plate_base == 384:
        plate_nrows, plate_ncols = 16, 24
    elif plate_base == 96:
        plate_nrows, plate_ncols = 8, 12
    else:
        raise ValueError(f'Plate base {plate_base} unknown')

    heatmap_data_list = []
    heatmap_names = []
    for plate, sub_df in data.groupby(groupby):
        # check if plate base are duplicated
        duplicated = sub_df[[row_base, col_base]].duplicated().sum() != 0
        if duplicated:
            if aggregation_func is None:
                raise ValueError('Row after groupby is not unique, aggregation_func can not be None')
            heatmap_data = sub_df.groupby([row_base, col_base])[value_col] \
                .apply(aggregation_func).unstack()
        else:
            heatmap_data = sub_df.set_index([row_base, col_base])[value_col] \
                .unstack()
        # reindex to make sure heatmap data in the shape of plate
        heatmap_data.index = range(heatmap_data.shape[0])
        heatmap_data.columns = range(heatmap_data.shape[1])
        heatmap_data = heatmap_data.reindex(index=list(range(plate_nrows)), columns=list(range(plate_ncols)))
        heatmap_data_list.append(heatmap_data)
        if isinstance(plate, str):
            heatmap_names.append(plate)
        else:
            heatmap_names.append('\n'.join(plate))

    nrows = round(len(heatmap_data_list) / ncols)
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1] * nrows),
                             ncols=ncols,
                             nrows=nrows)
    fig.suptitle(f'{value_col} on {len(heatmap_data_list)}*{plate_base} plates \n Color Range [{vmin}, {vmax}]',
                 fontsize=16)

    cmap = plt.cm.viridis
    cmap.set_under(color='#440154')
    cmap.set_over(color='#FDE725')
    cmap.set_bad(color='#FFFFFF')
    if heatmap_kws is None:
        heatmap_kws = {}
    for heatmap_data, heatmap_name, ax in zip(heatmap_data_list, heatmap_names, np.ravel(axes)):
        sns.heatmap(heatmap_data, vmin=vmin, vmax=vmax,
                    cmap=cmap, ax=ax, **heatmap_kws)
        ax.set(title=heatmap_name, ylim=(-0.5, plate_nrows + 0.5))
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    return fig, axes


def simple_violin(data, x, y, rotate_x=True):
    fig, ax = plt.subplots(figsize=(5, 3))
    sns.violinplot(x=x, y=y, data=data, ax=ax)
    if rotate_x:
        ax.tick_params(axis='x', rotation=90)
    return fig, ax


def cutoff_vs_cell_remain(data,
                          name='',
                          xlim_quantile=(0.001, 0.999),
                          ylim=None,
                          bins=100):
    xlim = tuple(np.quantile(data, xlim_quantile))
    x = np.linspace(xlim[0], xlim[1], 500)
    count_list = np.array([(data > i).sum() for i in x])
    original_total_data = data.size
    count_list = count_list / original_total_data * 100
    data = data[(data < xlim[1]) & (data > xlim[0])]

    fig, ax1 = plt.subplots()
    ax1 = sns.histplot(data, bins=bins, kde=False, ax=ax1)
    ax1.set_xlim(xlim)
    ax1.set_xlabel(name)
    if ylim is not None:
        ax1.set_ylim(*ylim)

    ax2 = ax1.twinx()
    ax2.plot(x,
             count_list,
             linewidth=2,
             linestyle='--',
             c='r')
    ax2.set_ylabel('% of Data Pass Filter', color='r')
    ax2.grid()
    return fig, (ax1, ax2)


def success_vs_fail(data, filter_col, filter_cutoff, x, y, ax):
    use_data = data.copy()
    use_data['filter'] = (data[filter_col] > filter_cutoff).apply(lambda i: 'Success' if i else 'Fail')
    sns.violinplot(x=x, y=y, data=use_data, hue='filter',
                   ax=ax)
    ax.tick_params(axis='x', rotation=90)
    return ax


def plot_dispersion(data, hue='gene_subset',
                    zlab='dispersion', data_quantile=(0.01, 0.99),
                    save_animate_path=None, fig_kws=None):
    from mpl_toolkits.mplot3d import Axes3D
    if Axes3D.__doc__:
        # touch the Axes3D to prevent ide remove it...
        pass

    @ticker.FuncFormatter
    def mean_formatter(x, pos):
        return f"{x:.1f}"

    _fig_kws = dict(figsize=(12, 4), dpi=160)
    if fig_kws is not None:
        _fig_kws.update(fig_kws)

    x = data['mean']
    y = data['cov']
    z = data[zlab]
    xlim = tuple(np.quantile(x, data_quantile))
    ylim = tuple(np.quantile(y, data_quantile))
    zlim = tuple(np.quantile(z, data_quantile))

    # directly apply lim on df
    _df = data[(x < xlim[1]) & (x > xlim[0]) &
               (y < ylim[1]) & (y > ylim[0]) &
               (z < zlim[1]) & (z > zlim[0])]
    color_dict = {True: 'steelblue',
                  False: 'lightgray'}
    color = _df[hue].map(color_dict).tolist()

    fig = plt.figure(**_fig_kws)
    if save_animate_path is None:
        axes = [fig.add_subplot(int(f'13{i}'), projection='3d') for i in range(1, 4)]
        view_inits = [(10, 10), (80, 45), (10, 80)]

        for i, (ax, view) in enumerate(zip(axes, view_inits)):
            ax.scatter(_df['mean'], _df['cov'], _df[zlab],
                       c=color, s=0.2, alpha=0.6)
            if i == 0:
                ax.set_xlabel('Mean', labelpad=10)
                ax.set_xticklabels([])
                ax.set_ylabel('Cov', labelpad=10)
                ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                ax.set_zlabel(zlab)
            elif i == 1:
                ax.set_xlabel('Mean', labelpad=10)
                ax.xaxis.set_major_formatter(mean_formatter)
                ax.set_ylabel('Cov', labelpad=10)
                ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                ax.set_zticklabels([])
            else:
                ax.set_xlabel('Mean', labelpad=10)
                ax.xaxis.set_major_formatter(mean_formatter)
                ax.set_yticklabels([])
                ax.set_ylabel('Cov', labelpad=10)
                ax.set_zlabel(zlab)
            ax.view_init(*view)
    else:
        axes = fig.add_subplot(111, projection='3d')
        axes.scatter(_df['mean'], _df['cov'], _df[zlab],
                     c=color, s=0.2, alpha=0.6)
        from matplotlib import animation

        def update(i):
            axes.view_init(10, i)

        ani = animation.FuncAnimation(fig, func=update,
                                      frames=100, interval=10, blit=False)
        ani.save(save_animate_path, writer='imagemagick')
    return fig, axes


def plot_hvf_selection():
    return
