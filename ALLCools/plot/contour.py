import numpy as np

from sklearn.neighbors import LocalOutlierFactor
from .utilities import zoom_min_max


def density_contour(ax, data, x, y, groupby=None, c='lightgray',
                    single_contour_pad=1, linewidth=1, palette=None):
    _data = data.copy()

    if groupby is not None:
        if isinstance(groupby, str):
            _data['groupby'] = data[groupby]
        else:
            _data['groupby'] = groupby
    else:
        _data['groupby'] = 'one group'

    _contour_kws = dict(linewidths=linewidth, levels=(-single_contour_pad, ), linestyles='dashed')
    _lof_kws = dict(n_neighbors=25, novelty=True, contamination='auto')

    xmin, ymin = _data[[x, y]].min()
    xmax, ymax = _data[[x, y]].max()
    xmin, xmax = zoom_min_max(xmin, xmax, 1.2)
    ymin, ymax = zoom_min_max(ymin, ymax, 1.2)

    for group, sub_data in _data[[x, y, 'groupby']].groupby('groupby'):
        xx, yy = np.meshgrid(np.linspace(xmin, xmax, 500), np.linspace(ymin, ymax, 500))
        clf = LocalOutlierFactor(**_lof_kws)
        clf.fit(sub_data.iloc[:, :2].values)
        z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
        z = z.reshape(xx.shape)
        if palette is None:
            _color = c
        else:
            _color = palette[group] if group in palette else c

        # plot contour line(s)
        ax.contour(xx, yy, z, colors=_color, **_contour_kws)
    return
