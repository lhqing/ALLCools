import seaborn as sns
from matplotlib.lines import Line2D

from .color import level_one_palette
from .contour import density_contour
from .text_anno_scatter import _text_anno_scatter
from .utilities import _make_tiny_axis_label, _density_based_sample, _extract_coords, zoom_ax


def categorical_scatter(
        data,
        ax,
        # coords
        coord_base='umap',
        x=None,
        y=None,
        # color
        hue=None,
        palette='auto',
        # text annotation
        text_anno=None,
        text_anno_kws=None,
        text_anno_palette=None,
        text_transform=None,
        dodge_text=False,
        dodge_kws=None,
        # legend
        show_legend=False,
        legend_kws=None,
        # size
        s=5,
        size=None,
        sizes: dict = None,
        size_norm=None,
        # other
        axis_format='tiny',
        max_points=5000,
        labelsize=4,
        linewidth=0,
        zoomxy=1.05,
        outline=None,
        outline_pad=3,
        outline_kws=None,
        scatter_kws=None,
):
    """
    Plot scatter plot with these options:
        - Color by a categorical variable, and generate legend of the variable if needed
        - Add text annotation using a categorical variable
        - Circle categories with outlines

    Parameters
    ----------
    data
        Dataframe that contains coordinates and categorical variables
    ax
        this function do not generate ax, must provide an ax
    coord_base
        coords name, if provided, will automatically search for x and y
    x
        x coord name
    y
        y coord name
    hue
        categorical col name or series for color hue
    palette
        palette of the hue, str or dict
    text_anno
        categorical col name or series for text annotation
    text_anno_kws
    text_anno_palette
    text_transform
    dodge_text
    dodge_kws
    show_legend
    legend_kws
    s
    size
    sizes
    size_norm
    axis_format
    max_points
    labelsize
    linewidth
    zoomxy
    outline
    outline_pad
    outline_kws
    scatter_kws
        kws dict pass to sns.scatterplot

    Returns
    -------

    """
    # add coords
    _data, x, y = _extract_coords(data, coord_base, x, y)
    # _data has 2 cols: "x" and "y"

    # down sample plot data if needed.
    if max_points is not None:
        if _data.shape[0] > max_points:
            _data = _density_based_sample(_data, seed=1, size=max_points,
                                          coords=['x', 'y'])

    # default scatter options
    _scatter_kws = {'linewidth': 0, 's': s, 'legend': None, 'palette': palette}
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # deal with color
    palette_dict = None
    if hue is not None:
        if isinstance(hue, str):
            _data['hue'] = data[hue]
        else:
            _data['hue'] = hue
        hue = 'hue'
        _data['hue'] = _data['hue'].astype('category')
        # deal with color palette
        palette = _scatter_kws['palette']
        if isinstance(palette, str) or isinstance(palette, list):
            palette_dict = level_one_palette(_data['hue'], order=None, palette=palette)
        elif isinstance(palette, dict):
            palette_dict = palette
        else:
            raise TypeError(f'Palette can only be str, list or dict, '
                            f'got {type(palette)}')
        _scatter_kws['palette'] = palette_dict

    # deal with size
    if size is not None:
        # discard s from _scatter_kws and use size in sns.scatterplot
        s = _scatter_kws.pop('s')
    sns.scatterplot(x='x', y='y', data=_data, ax=ax, hue=hue,
                    size=size, sizes=sizes, size_norm=size_norm,
                    **_scatter_kws)

    # deal with text annotation
    if text_anno:
        if isinstance(text_anno, str):
            _data['text_anno'] = data[text_anno]
        else:
            _data['text_anno'] = text_anno

        _text_anno_scatter(data=_data[['x', 'y', 'text_anno']],
                           ax=ax,
                           x='x',
                           y='y',
                           dodge_text=dodge_text,
                           dodge_kws=dodge_kws,
                           palette=text_anno_palette,
                           text_transform=text_transform,
                           anno_col='text_anno',
                           text_anno_kws=text_anno_kws,
                           labelsize=labelsize)

    # deal with outline
    if outline:
        if isinstance(outline, str):
            _data['outline'] = data[outline]
        else:
            _data['outline'] = outline
        _outline_kws = {'linewidth': linewidth,
                        'palette': None,
                        'c': 'lightgray',
                        'single_contour_pad': outline_pad}
        if outline_kws is not None:
            _outline_kws.update(outline_kws)
        density_contour(ax=ax, data=_data, x='x', y='y', groupby='outline', **_outline_kws)

    # clean axis
    if axis_format == 'tiny':
        _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
    elif (axis_format == 'empty') or (axis_format is None):
        sns.despine(ax=ax, left=True, bottom=True)
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    else:
        pass

    # deal with legend
    if show_legend and (hue is not None):
        n_hue = len(palette_dict)
        _legend_kws = dict(ncol=(1 if n_hue <= 20 else 2 if n_hue <= 40 else 3),
                           fontsize=labelsize,
                           bbox_to_anchor=(1.05, 1),
                           loc='upper left',
                           borderaxespad=0.)
        if legend_kws is not None:
            _legend_kws.update(legend_kws)

        handles = []
        labels = []
        exist_hues = _data['hue'].unique()
        for hue_name, color in palette_dict.items():
            if hue_name not in exist_hues:
                # skip hue_name that do not appear in the plot
                continue
            handle = Line2D([0], [0], marker='o', color='w',
                            markerfacecolor=color, markersize=_legend_kws['fontsize'])
            handles.append(handle)
            labels.append(hue_name)
        _legend_kws['handles'] = handles
        _legend_kws['labels'] = labels
        ax.legend(**_legend_kws)

    if zoomxy is not None:
        zoom_ax(ax, zoomxy)

    return ax, _data
