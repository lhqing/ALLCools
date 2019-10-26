import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

from .color import level_one_palette
from .text_anno_scatter import text_anno_scatter
from .utilities import _make_tiny_axis_label, _density_based_sample


def categorical_scatter(data, ax, coord_base='umap', x=None, y=None,
                        scatter_kws=None,  # about basic scatter
                        hue=None, palette='tab10',  # about hue
                        text_anno=None, dodge_text=False, dodge_kws=None,  # about text anno
                        text_anno_kws=None, text_anno_palette=None,  # about text anno
                        show_legend=False, legend_kws=None,  # about legend
                        axis_format='tiny', max_points=5000, labelsize=4):  # other adjustment
    data = data.copy()
    # down sample plot data if needed.
    if max_points is not None:
        if data.shape[0] > max_points:
            data = _density_based_sample(data, seed=1, size=max_points,
                                         coords=[f'{coord_base}_0',
                                                 f'{coord_base}_1'])

    # default scatter options
    _scatter_kws = {'linewidth': 0, 's': 5, 'legend': None, 'palette': palette}
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # add coords as x and y
    if coord_base is not None:
        x = f'{coord_base}_0'
        y = f'{coord_base}_1'
        if (x not in data.columns) or (y not in data.columns):
            raise KeyError(f'coord_base {coord_base} is provided, '
                           f'but {coord_base}_0 or {coord_base}_1 not found in columns')
    else:
        if (x is None) or (y is None):
            raise ValueError('Either provide coord_base, or provide both x and y.')
    _data = pd.DataFrame({'x': data[x],
                          'y': data[y]})

    palette_dict = None
    if hue is not None:
        if isinstance(hue, str):
            _data[hue] = data[hue].astype('category')
        else:
            _data['hue'] = hue.astype('category')
            hue = 'hue'
        # deal with color palette
        palette = _scatter_kws['palette']
        if isinstance(palette, str) or isinstance(palette, list):
            palette_dict = level_one_palette(_data[hue], order=None, palette=palette)
        elif isinstance(palette, dict):
            palette_dict = palette
        else:
            raise TypeError(f'Palette can only be str, list or dict, '
                            f'got {type(palette)}')
        _scatter_kws['palette'] = palette_dict

    sns.scatterplot(x='x', y='y', hue=hue,
                    data=_data, ax=ax, **_scatter_kws)

    # clean axis
    if axis_format == 'tiny':
        _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
    elif (axis_format == 'empty') or (axis_format is None):
        sns.despine(ax=ax, left=True, bottom=True)
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    else:
        pass

    if text_anno:
        if isinstance(text_anno, str):
            _data['text_anno'] = data[text_anno]
        else:
            _data['text_anno'] = text_anno

        text_anno_scatter(data=_data[['x', 'y', 'text_anno']],
                          ax=ax,
                          x='x',
                          y='y',
                          dodge_text=dodge_text,
                          dodge_kws=dodge_kws,
                          palette=text_anno_palette,
                          anno_col=text_anno,
                          text_anno_kws=text_anno_kws)

    if show_legend and (hue is not None):
        n_hue = len(palette_dict)
        _legend_kws = dict(ncol=(1 if n_hue <= 14 else 2 if n_hue <= 30 else 3), fontsize=labelsize)
        if legend_kws is not None:
            _legend_kws.update(legend_kws)

        handles = []
        labels = []
        exist_hues = _data[hue].unique()
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
    return ax, _data
