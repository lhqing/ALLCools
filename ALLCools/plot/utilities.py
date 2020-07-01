from decimal import Decimal

import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.neighbors import LocalOutlierFactor


def _density_based_sample(data: pd.DataFrame, coords: list, portion=None, size=None, seed=None):
    """down sample data based on density, to prevent overplot in dense region and decrease plotting time"""
    clf = LocalOutlierFactor(n_neighbors=20, algorithm='auto',
                             leaf_size=30, metric='minkowski',
                             p=2, metric_params=None, contamination=0.1)

    # coords should already exist in data, get them by column names list
    data_coords = data[coords]
    clf.fit(data_coords)
    # original score is negative, the larger the denser
    density_score = clf.negative_outlier_factor_
    delta = density_score.max() - density_score.min()
    # density score to probability: the denser the less probability to be picked up
    probability_score = 1 - (density_score - density_score.min()) / delta
    probability_score = np.sqrt(probability_score)
    probability_score = probability_score / probability_score.sum()

    if size is not None:
        pass
    elif portion is not None:
        size = int(data_coords.index.size * portion)
    else:
        raise ValueError('Either portion or size should be provided.')
    if seed is not None:
        np.random.seed(seed)
    selected_cell_index = np.random.choice(data_coords.index,
                                           size=size,
                                           replace=False,
                                           p=probability_score)  # choice data based on density weights

    # return the down sampled data
    return data.reindex(selected_cell_index)


def _translate_coord_name(coord_name):
    return coord_name.upper().replace('_', ' ')


def _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5):
    """this function assume coord is [0, 1]"""
    # clean ax axises
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    sns.despine(ax=ax, left=True, bottom=True)

    _arrow_kws = dict(width=0.003, linewidth=0, color='black')
    if arrow_kws is not None:
        _arrow_kws.update(arrow_kws)

    ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws,
             transform=ax.transAxes)
    ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws,
             transform=ax.transAxes)
    ax.text(0.06, 0.03, _translate_coord_name(x),
            fontdict=dict(fontsize=fontsize,
                          horizontalalignment='left',
                          verticalalignment='center'),
            transform=ax.transAxes)
    ax.text(0.03, 0.06, _translate_coord_name(y),
            fontdict=dict(fontsize=fontsize,
                          rotation=90, rotation_mode='anchor',
                          horizontalalignment='left',
                          verticalalignment='center'),
            transform=ax.transAxes)
    return


def _extract_coords(data, coord_base, x, y):
    if (x is not None) and (y is not None):
        pass
    else:
        x = f'{coord_base}_0'
        y = f'{coord_base}_1'
    if (x not in data.columns) or (y not in data.columns):
        raise KeyError(f'{x} or {y} not found in columns.')

    _data = pd.DataFrame({'x': data[x],
                          'y': data[y]})
    return _data, x, y


def zoom_min_max(vmin, vmax, scale):
    width = vmax - vmin
    width_zoomed = width * scale
    delta_value = (width_zoomed - width) / 2
    return vmin - delta_value, vmax + delta_value


def zoom_ax(ax, zoom_scale, on='both'):
    on = on.lower()

    xlim = ax.get_xlim()
    xlim_zoomed = zoom_min_max(*xlim, zoom_scale)

    ylim = ax.get_ylim()
    ylim_zoomed = zoom_min_max(*ylim, zoom_scale)

    if (on == 'both') or ('x' in on):
        ax.set_xlim(xlim_zoomed)
    if (on == 'both') or ('y' in on):
        ax.set_ylim(ylim_zoomed)


def smart_number_format(x, pos=None):
    if (x > 0.01) and (x < 1):
        return f'{x:.2f}'.rstrip('0')
    elif (x >= 1) and (x < 100):
        return f'{int(x)}'
    else:
        t = f"{Decimal(x):.2E}"
        if t == '0.00E+2':
            return '0'
        else:
            return t


def add_ax_box(ax, expend=0, **patch_kws):
    import matplotlib.patches as patches
    _patch_kws = dict(linewidth=1,
                      edgecolor='k',
                      facecolor='none')
    _patch_kws.update(patch_kws)

    rect = patches.Rectangle((0 - expend, 0 - expend),
                             1 + expend,
                             1 + expend,
                             transform=ax.transAxes,
                             **_patch_kws)

    # Add the patch to the Axes
    ax.add_patch(rect)
    return ax


def tight_hue_range(hue_data, portion):
    """Automatic select a SMALLEST data range that covers [portion] of the data"""
    hue_data = hue_data[np.isfinite(hue_data)]
    hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
    min_window_right = hue_quantiles.rolling(window=int(portion * 100)) \
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
    return vmin, vmax
