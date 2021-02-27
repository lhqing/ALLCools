import numpy as np
import seaborn as sns
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import copy
from .color import plot_colorbar
from .contour import density_contour
from .text_anno_scatter import _text_anno_scatter
from .utilities import _density_based_sample, _extract_coords, _make_tiny_axis_label, zoom_ax


def tight_hue_range(hue_data, portion):
    """Automatic select a SMALLEST data range that covers [portion] of the data"""
    hue_data = hue_data[np.isfinite(hue_data)]
    hue_quantile = hue_data.quantile(q=np.arange(0, 1, 0.01))
    min_window_right = hue_quantile.rolling(window=int(portion * 100)) \
        .apply(lambda i: i.max() - i.min(), raw=True) \
        .idxmin()
    min_window_left = max(0, min_window_right - portion)
    vmin, vmax = tuple(hue_data.quantile(q=[min_window_left,
                                            min_window_right]))
    if np.isfinite(vmin):
        vmin = max(hue_data.min(), vmin)
    else:
        vmin = hue_data.min()
    if np.isfinite(vmax):
        vmax = min(hue_data.max(), vmax)
    else:
        vmax = hue_data.max()

    if vmin == vmax:
        return hue_data.min, hue_data.max
    return vmin, vmax


def continuous_scatter(
        data,
        ax,
        coord_base='umap',
        x=None,
        y=None,
        scatter_kws=None,
        hue=None,
        hue_norm=None,
        hue_portion=0.95,
        cmap='viridis',
        colorbar=True,
        colorbar_label_kws=None,
        size=None,
        size_norm=None,
        size_portion=0.95,
        sizes=None,
        sizebar=True,
        text_anno=None,
        dodge_text=False,
        dodge_kws=None,
        text_anno_kws=None,
        text_anno_palette=None,
        text_transform=None,
        axis_format='tiny',
        max_points=5000,
        s=5,
        labelsize=4,
        linewidth=.5,
        cax=None,
        zoomxy=1.05,
        outline=None,
        outline_kws=None,
        outline_pad=2
):
    # add coords
    _data, x, y = _extract_coords(data, coord_base, x, y)
    # _data has 2 cols: "x" and "y"

    # down sample plot data if needed.
    if max_points is not None:
        if _data.shape[0] > max_points:
            _data = _density_based_sample(_data, seed=1, size=max_points,
                                          coords=['x', 'y'])

    # default scatter options
    _scatter_kws = {'linewidth': 0, 's': s, 'legend': None}
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # deal with color
    if hue is not None:
        if isinstance(hue, str):
            _data['hue'] = data[hue]
            colorbar_label = hue
        else:
            _data['hue'] = hue
            colorbar_label = hue.name
        _data['hue'] = _data['hue'].astype(float)
        hue = 'hue'

        if hue_norm is None:
            # get the smallest range that include "hue_portion" of data
            hue_norm = tight_hue_range(_data['hue'], hue_portion)
        if isinstance(cmap, str):
            # from here, cmap become colormap object
            cmap = copy.copy(get_cmap(cmap))
            cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5))
        else:
            if not isinstance(cmap, ScalarMappable):
                raise TypeError(f'cmap can only be str or ScalarMappable, got {type(cmap)}')
        # cnorm is the normalizer for color
        cnorm = Normalize(vmin=hue_norm[0],
                          vmax=hue_norm[1])
    else:
        hue_norm = None
        cnorm = None
        colorbar_label = ''

    # deal with size
    if size is not None:
        if isinstance(size, str):
            _data['size'] = data[size].astype(float)
        else:
            _data['size'] = size.astype(float)
        size = 'size'

        if size_norm is None:
            # get the smallest range that include "size_portion" of data
            size_norm = tight_hue_range(_data['size'], size_portion)

            # snorm is the normalizer for size
            size_norm = Normalize(vmin=size_norm[0],
                                  vmax=size_norm[1])
        # replace s with sizes
        s = _scatter_kws.pop('s')
        if sizes is None:
            sizes = (min(s, 1), s)
    else:
        size_norm = None
        sizes = None

    sns.scatterplot(x='x', y='y', data=_data,
                    hue=hue, palette=cmap, hue_norm=cnorm,
                    size=size, sizes=sizes, size_norm=size_norm,
                    ax=ax, **_scatter_kws)

    if text_anno is not None:
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

    return_axes = [ax]

    # make color bar
    if colorbar and (hue is not None):
        _colorbar_label_kws = dict(fontsize=labelsize, label=hue, labelpad=10, rotation=270)
        if colorbar_label_kws is not None:
            _colorbar_label_kws.update(colorbar_label_kws)

        # small ax for colorbar
        if cax is None:
            cax = inset_axes(ax, width="3%", height="25%",
                             loc='lower right', borderpad=0)
        cax = plot_colorbar(cax=cax, cmap=cmap, cnorm=cnorm, hue_norm=hue_norm,
                            label=colorbar_label, orientation='vertical', labelsize=labelsize, linewidth=0.5)
        return_axes.append(cax)

    # make size bar
    if sizebar and (size is not None):
        # TODO plot dot size bar
        pass

    if zoomxy is not None:
        zoom_ax(ax, zoomxy)

    return tuple(return_axes), _data
