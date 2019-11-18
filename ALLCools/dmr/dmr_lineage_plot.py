import seaborn as sns
from matplotlib.cm import get_cmap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram


def straight_branch(ax, a, b, plot_kws):
    """Draw link line between a and b"""
    a_x, ay = a
    bx, by = b
    branch_x = [a_x, bx, bx]
    branch_y = [ay, ay, by]
    if plot_kws is None:
        plot_kws = {}
    return ax.plot(branch_x, branch_y, **plot_kws)


def plot_dendrogram(linkage_df,
                    labels_list,
                    dendro_kws=None,
                    ax=None,
                    branch_type='straight',
                    plot_node_id=False,
                    plot_kws=None,
                    node_hue=None,
                    node_hue_norm=None,
                    palette='viridis',  # shared by both line and node hue
                    node_size=None,
                    node_size_norm=None,
                    line_hue=None,
                    line_hue_norm=None,
                    sizes=None,
                    size=30,
                    linewidth=1,
                    color=None):
    if plot_kws is None:
        plot_kws = {}

    _dendro_kws = dict(no_plot=True)
    if dendro_kws is not None:
        _dendro_kws.update(dendro_kws)
    # all we need is the leaves order from dendrogram,
    # bellow we recalculate the node position to match the node id,
    # so we can control each node
    dendro = dendrogram(linkage_df, labels=labels_list, **_dendro_kws)
    n_leaves = len(dendro['leaves'])

    node_pos = {}  # all node including singleton and non-singleton
    for leaf_x, leaf in enumerate(dendro['leaves']):
        node_pos[int(leaf)] = (leaf_x, 0)

    direct_link_map = {}  # node linkage, keys only contain non-singleton
    for i, (left, right, height, _) in linkage_df.iterrows():
        node_id = int(i + linkage_df.shape[0] + 1)
        left = int(left)
        right = int(right)
        node_x = (node_pos[left][0] + node_pos[right][0]) / 2
        node_pos[node_id] = [node_x, height]
        direct_link_map[node_id] = [int(left), int(right)]

    if branch_type == 'straight':
        branch_plot_function = straight_branch
    else:
        raise

    # node colors
    nan_color = '#D3D3D3' if color is None else color
    if node_hue is not None:
        if node_hue_norm is None:
            values = node_hue.values()
            _hue_norm = Normalize(vmin=min(values),
                                  vmax=max(values))
        else:
            _hue_norm = Normalize(vmin=min(node_hue_norm),
                                  vmax=max(node_hue_norm))
        _cmap = get_cmap(palette)

        def node_cmap(v):
            return (_cmap(_hue_norm(v)),)
    else:
        node_hue = {}

        def node_cmap(_):
            raise
    node_colors = {node: node_cmap(node_hue[node]) if (node in node_hue) else nan_color
                   for node in node_pos.keys()}

    # node sizes
    nan_size = size
    if node_size is not None:
        if sizes is None:
            nan_size = 10
            sizes = (10, 80)
        else:
            nan_size = sizes[0]

        if node_size_norm is None:
            values = node_size.values()
            _size_norm = Normalize(vmin=min(values),
                                   vmax=max(values))
        else:
            _size_norm = Normalize(vmin=min(node_size_norm),
                                   vmax=max(node_size_norm))

        def node_smap(v):
            v_norm = _size_norm(v)
            v_norm = min(1, max(0, v_norm)) # limit norm value to [0, 1]
            s = v_norm * (sizes[1] - sizes[0]) + sizes[0]
            return s
    else:
        node_size = {}

        def node_smap(_):
            raise
    node_sizes = {node: node_smap(node_size[node]) if (node in node_size) else nan_size
                  for node in node_pos.keys()}

    for node_id, (node_x, node_y) in node_pos.items():
        ax.scatter(node_x, node_y, s=node_sizes[node_id], c=node_colors[node_id], zorder=3)

    # line color
    nan_color = '#D3D3D3' if color is None else color
    if line_hue is not None:
        if line_hue_norm is None:
            values = line_hue.values()
            _hue_norm = Normalize(vmin=min(values),
                                  vmax=max(values))
        else:
            _hue_norm = Normalize(vmin=min(line_hue_norm),
                                  vmax=max(line_hue_norm))
        _cmap = get_cmap(palette)

        def line_cmap(v):
            return _cmap(_hue_norm(v))
    else:
        line_hue = {}

        def line_cmap(_):
            raise
    line_colors = {node: line_cmap(line_hue[node]) if (node in line_hue) else nan_color
                   for node in node_pos.keys()}

    ymax = 0
    for node_id, (node_x, node_y) in node_pos.items():
        ymax = max(ymax, node_y)
        # plot node id text
        if plot_node_id:
            if node_id >= n_leaves:
                ax.text(node_x, node_y, node_id,
                        fontsize=4 if 'fontsize' not in plot_kws else plot_kws['fontsize'],
                        ha='center', va='center', c='k')
            else:
                ax.text(node_x, -0.01, node_id,
                        fontsize=4 if 'fontsize' not in plot_kws else plot_kws['fontsize'],
                        ha='center', va='center', c='k')

        # plot branch
        # only non-singleton node has branch:
        if node_id in direct_link_map:
            # get child
            left_child, right_child = direct_link_map[node_id]

            # plot left branch
            branch_plot_function(ax, (node_x, node_y),
                                 node_pos[left_child],
                                 plot_kws=dict(c=line_colors[left_child],
                                               linewidth=linewidth))

            # plot right branch
            branch_plot_function(ax, (node_x, node_y),
                                 node_pos[right_child],
                                 plot_kws=dict(c=line_colors[right_child],
                                               linewidth=linewidth))

    ax.set_ylim(0 - ymax * 0.05, ymax * 1.05)
    return node_pos, dendro


def plot_parsimony_data(data, ax, hue_norm, palette):
    sns.barplot(data=data, x='index', y='tree',
                label='tree', ax=ax,
                color='steelblue', alpha=0.3)
    sns.scatterplot(data=data, x='index', y='raw',
                    label='raw', ax=ax, color='lightgray',
                    hue='tree', hue_norm=hue_norm, palette=palette,
                    alpha=0.7, linewidth=1, legend=None)
    return


def plot_colorbar(cax, cmap, cnorm, hue_norm, linewidth=0.5):
    if isinstance(cmap, str):
        cmap = get_cmap(cmap)

    colorbar = ColorbarBase(cax, cmap=cmap, norm=cnorm,
                            orientation='vertical', extend='both')
    colorbar_ticks = [hue_norm[0], sum(hue_norm) / 2, hue_norm[1]]
    # TODO automatic ticklabel format, auto sci-format, float trim etc
    colorbar_ticklabels = list(map(lambda i: f'{i:.1f}', colorbar_ticks))
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels(colorbar_ticklabels)
    colorbar.outline.set_linewidth(linewidth)
    colorbar.ax.tick_params(size=1, pad=1, width=linewidth, length=0.3)
    return cax
