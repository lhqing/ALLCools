import seaborn as sns
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram as _dendrogram
import matplotlib as mpl


def straight_branch(ax, a, b, plot_kws):
    """Draw link line between a and b"""
    a_x, ay = a
    bx, by = b
    branch_x = [a_x, bx, bx]
    branch_y = [ay, ay, by]
    if plot_kws is None:
        plot_kws = {}
    return ax.plot(branch_x, branch_y, **plot_kws)


def plot_dendrogram(
        linkage_df,
        ax,
        dendro=None,
        labels=None,
        dendro_kws=None,
        plot_node_id=False,
        plot_non_singleton=True,
        plot_kws=None,
        node_hue=None,
        node_hue_norm=None,
        node_hue_cbar=True,
        node_hue_cbar_frac=0.1,
        node_palette='viridis',  # shared by both line and node hue
        node_size=None,
        node_size_norm=None,
        node_sizes=None,
        line_hue=None,
        line_hue_norm=None,
        line_palette='gray_r',
        linewidth=1.5,
        edge_color='gray',
        marker_size=60,
        marker_color='lightblue'):
    """

    Parameters
    ----------
    linkage_df
    dendro
    labels
    dendro_kws
    ax
    plot_node_id
    plot_non_singleton
    plot_kws
    node_hue
    node_hue_norm
    node_hue_cbar
    node_hue_cbar_frac
    node_palette
    node_size
    node_size_norm
    node_sizes
    line_hue
    line_hue_norm
    line_palette
    linewidth
    edge_color
    marker_size
    marker_color

    Returns
    -------

    """
    if plot_kws is None:
        plot_kws = {}

    if dendro is None:
        if labels is None or linkage_df is None:
            raise ValueError(
                'linkage_df and labels must be provided to calculate dendrogram.')
        print('Computing dendrogram')
        _dendro_kws = dict(no_plot=True)
        if dendro_kws is not None:
            _dendro_kws.update(dendro_kws)
        # all we need is the leaves order from dendrogram,
        # bellow we recalculate the node position to match the node id,
        # so we can control each node
        dendro = _dendrogram(linkage_df, labels=labels, **_dendro_kws)
    else:
        labels = dendro['ivl']
    n_leaves = len(dendro['leaves'])

    node_pos = {}  # all node including singleton and non-singleton
    direct_link_map = {}  # node linkage, keys only contain non-singleton
    for leaf_x, leaf in enumerate(dendro['leaves']):
        # add singleton positions first
        node_pos[int(leaf)] = (leaf_x, 0)
    for i, (idx, (left, right, height, _)) in enumerate(linkage_df.iterrows()):
        node_id = int(i + linkage_df.shape[0] + 1)
        left = int(left)
        right = int(right)
        node_x = (node_pos[left][0] + node_pos[right][0]) / 2
        node_pos[node_id] = [node_x, height]
        direct_link_map[node_id] = [int(left), int(right)]

    # ------------------ Plot nodes ------------------
    # node colors
    nan_color = '#D3D3D3' if marker_color is None else marker_color
    if node_hue is not None:
        if node_hue_norm is None:
            values = node_hue.values
            _hue_norm = Normalize(vmin=min(values), vmax=max(values))
        else:
            _hue_norm = Normalize(vmin=min(node_hue_norm),
                                  vmax=max(node_hue_norm))
        _cmap = get_cmap(node_palette)

        def node_cmap(v):
            return (_cmap(_hue_norm(v)),)

        if node_hue_cbar:
            ax.figure.colorbar(mpl.cm.ScalarMappable(norm=_hue_norm, cmap=_cmap),
                               ax=ax.figure.axes,
                               shrink=0.6,
                               fraction=node_hue_cbar_frac,
                               label='Node Color')
    else:
        node_hue = {}

        def node_cmap(_):
            return nan_color

    node_colors = {
        node: node_cmap(node_hue[node]) if (node in node_hue) else nan_color
        for node in node_pos.keys()
    }
    # node sizes
    nan_size = marker_size
    if node_size is not None:
        if node_sizes is None:
            node_sizes = (marker_size, marker_size * 2)
        if node_size_norm is None:
            values = node_size.values
            _size_norm = Normalize(vmin=min(values), vmax=max(values))
        else:
            _size_norm = Normalize(vmin=min(node_size_norm),
                                   vmax=max(node_size_norm))

        def node_smap(v):
            v_norm = _size_norm(v)
            v_norm = min(1, max(0, v_norm))  # limit norm value to [0, 1]
            s = v_norm * (node_sizes[1] - node_sizes[0]) + node_sizes[0]
            return s
    else:
        node_size = {}

        def node_smap(_):
            return nan_size

    node_sizes = {
        node: node_smap(node_size[node]) if (node in node_size) else nan_size
        for node in node_pos.keys()
    }
    # plot nodes
    for node_id, (node_x, node_y) in node_pos.items():
        if (node_id > len(dendro['leaves']) - 1) and not plot_non_singleton:
            break
        ax.scatter(node_x,
                   node_y,
                   s=node_sizes[node_id],
                   c=node_colors[node_id],
                   clip_on=False,
                   zorder=3)

    # ------------------ Plot edges and node id ------------------
    # line color
    nan_color = '#D3D3D3' if edge_color is None else edge_color
    if line_hue is not None:
        if line_hue_norm is None:
            values = line_hue.values
            _hue_norm = Normalize(vmin=min(values), vmax=max(values))
        else:
            _hue_norm = Normalize(vmin=min(line_hue_norm),
                                  vmax=max(line_hue_norm))
        _cmap = get_cmap(line_palette)

        def line_cmap(v):
            return _cmap(_hue_norm(v))
    else:
        line_hue = {}

        def line_cmap(_):
            return nan_color

    line_colors = {
        node: line_cmap(line_hue[node]) if (node in line_hue) else nan_color
        for node in node_pos.keys()
    }

    ymax = 0
    for node_id, (node_x, node_y) in node_pos.items():
        ymax = max(ymax, node_y)
        # plot node id text
        if plot_node_id:
            if node_id >= n_leaves:
                ax.text(node_x,
                        node_y,
                        node_id,
                        fontsize=6 if 'fontsize' not in plot_kws else
                        plot_kws['fontsize'],
                        ha='center',
                        va='center',
                        c='k')
            else:
                ax.text(node_x,
                        -0.01,
                        node_id,
                        fontsize=6 if 'fontsize' not in plot_kws else
                        plot_kws['fontsize'],
                        ha='center',
                        va='center',
                        c='k')

        # plot branch
        # only non-singleton node has branch:
        if node_id in direct_link_map:
            # get child
            left_child, right_child = direct_link_map[node_id]
            # plot left branch
            straight_branch(ax, (node_x, node_y),
                            node_pos[left_child],
                            plot_kws=dict(c=line_colors[left_child],
                                          linewidth=linewidth,
                                          clip_on=False))
            # plot right branch
            straight_branch(ax, (node_x, node_y),
                            node_pos[right_child],
                            plot_kws=dict(c=line_colors[right_child],
                                          linewidth=linewidth,
                                          clip_on=False))

    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set(xticks=range(len(labels)), xticklabels=dendro['ivl'])
    ax.xaxis.set_tick_params(rotation=90)
    sns.despine(ax=ax, bottom=True, offset=15)
    return dendro, node_pos
