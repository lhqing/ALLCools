from decimal import Decimal

import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from sklearn.neighbors import LocalOutlierFactor


def _density_based_sample(data: pd.DataFrame, coords: list, portion=None, size=None, seed=None):
    """Down sample data based on density, to prevent overplot in dense region and decrease plotting time."""
    clf = LocalOutlierFactor(
        n_neighbors=20,
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        metric_params=None,
        contamination=0.1,
    )

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
        raise ValueError("Either portion or size should be provided.")
    if seed is not None:
        np.random.seed(seed)
    selected_cell_index = np.random.choice(
        data_coords.index, size=size, replace=False, p=probability_score
    )  # choice data based on density weights

    # return the down sampled data
    return data.reindex(selected_cell_index)


def _translate_coord_name(coord_name):
    return coord_name.upper().replace("_", " ")


def _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5):
    # This function assume coord is [0, 1].

    # clean ax axises
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    sns.despine(ax=ax, left=True, bottom=True)

    _arrow_kws = {"width": 0.003, "linewidth": 0, "color": "black"}
    if arrow_kws is not None:
        _arrow_kws.update(arrow_kws)

    ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws, transform=ax.transAxes)
    ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws, transform=ax.transAxes)
    ax.text(
        0.06,
        0.03,
        _translate_coord_name(x),
        fontdict={"fontsize": fontsize, "horizontalalignment": "left", "verticalalignment": "center"},
        transform=ax.transAxes,
    )
    ax.text(
        0.03,
        0.06,
        _translate_coord_name(y),
        fontdict={
            "fontsize": fontsize,
            "rotation": 90,
            "rotation_mode": "anchor",
            "horizontalalignment": "left",
            "verticalalignment": "center",
        },
        transform=ax.transAxes,
    )
    return


def _extract_coords(data, coord_base, x, y):
    if (x is not None) and (y is not None):
        pass
    else:
        x = f"{coord_base}_0"
        y = f"{coord_base}_1"

    if isinstance(data, anndata.AnnData):
        adata = data
        _data = pd.DataFrame(
            {
                "x": adata.obsm[f"X_{coord_base}"][:, 0],
                "y": adata.obsm[f"X_{coord_base}"][:, 1],
            },
            index=adata.obs_names,
        )
    elif isinstance(data, xr.Dataset):
        ds = data
        if coord_base not in ds.dims:
            raise KeyError(f"xr.Dataset do not contain {coord_base} dim")
        data_var = {i for i in ds.data_vars.keys() if i.startswith(coord_base)}.pop()
        _data = pd.DataFrame(
            {
                "x": ds[data_var].sel({coord_base: f"{coord_base}_0"}).to_pandas(),
                "y": ds[data_var].sel({coord_base: f"{coord_base}_1"}).to_pandas(),
            }
        )
    else:
        if (x not in data.columns) or (y not in data.columns):
            raise KeyError(f"{x} or {y} not found in columns.")
        _data = pd.DataFrame({"x": data[x], "y": data[y]})
    return _data, x, y


def zoom_min_max(vmin, vmax, scale):
    """Zoom min and max value."""
    width = vmax - vmin
    width_zoomed = width * scale
    delta_value = (width_zoomed - width) / 2
    return vmin - delta_value, vmax + delta_value


def zoom_ax(ax, zoom_scale, on="both"):
    """Zoom ax on both x and y-axis."""
    on = on.lower()

    xlim = ax.get_xlim()
    xlim_zoomed = zoom_min_max(*xlim, zoom_scale)

    ylim = ax.get_ylim()
    ylim_zoomed = zoom_min_max(*ylim, zoom_scale)

    if (on == "both") or ("x" in on):
        ax.set_xlim(xlim_zoomed)
    if (on == "both") or ("y" in on):
        ax.set_ylim(ylim_zoomed)


def smart_number_format(x, pos=None):
    """Format numbers automatically."""
    if (x > 0.01) and (x < 1):
        return f"{x:.2f}".rstrip("0")
    elif (x >= 1) and (x < 100):
        return f"{int(x)}"
    else:
        t = f"{Decimal(x):.2E}"
        if t == "0.00E+2":
            return "0"
        else:
            return t


def add_ax_box(ax, expend=0, **patch_kws):
    """Add a box around the ax."""
    from matplotlib import patches

    _patch_kws = {"linewidth": 1, "edgecolor": "k", "facecolor": "none"}
    _patch_kws.update(patch_kws)

    rect = patches.Rectangle(
        (0 - expend, 0 - expend),
        1 + expend,
        1 + expend,
        transform=ax.transAxes,
        **_patch_kws,
    )

    # Add the patch to the Axes
    ax.add_patch(rect)
    return ax


def tight_hue_range(hue_data, portion):
    """Automatic select a SMALLEST data range that covers [portion] of the data."""
    hue_data = hue_data[np.isfinite(hue_data)]
    hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
    min_window_right = (
        hue_quantiles.rolling(window=int(portion * 100)).apply(lambda i: i.max() - i.min(), raw=True).idxmin()
    )
    min_window_left = max(0, min_window_right - portion)
    vmin, vmax = tuple(hue_data.quantile(q=[min_window_left, min_window_right]))
    if np.isfinite(vmin):
        vmin = max(hue_data.min(), vmin)
    else:
        vmin = hue_data.min()
    if np.isfinite(vmax):
        vmax = min(hue_data.max(), vmax)
    else:
        vmax = hue_data.max()
    return vmin, vmax


def _take_data_series(data, k):
    if isinstance(data, (xr.Dataset, xr.DataArray)):
        _value = data[k].to_pandas()
    elif isinstance(data, anndata.AnnData):
        _value = data.obs[k].copy()
    else:
        _value = data[k].copy()
    return _value


def _auto_size(ax, n_dots):
    """Auto determine dot size based on ax size and n dots"""
    bbox = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted())
    scale = bbox.width * bbox.height / 14.6  # 14.6 is a 5*5 fig I used to estimate
    n = n_dots / scale  # larger figure means data look sparser
    if n < 500:
        s = 14 - n / 100
    elif n < 1500:
        s = 7
    elif n < 3000:
        s = 5
    elif n < 8000:
        s = 3
    elif n < 15000:
        s = 2
    elif n < 30000:
        s = 1.5
    elif n < 50000:
        s = 1
    elif n < 80000:
        s = 0.8
    elif n < 150000:
        s = 0.6
    elif n < 300000:
        s = 0.5
    elif n < 500000:
        s = 0.4
    elif n < 800000:
        s = 0.3
    elif n < 1000000:
        s = 0.2
    elif n < 2000000:
        s = 0.1
    elif n < 3000000:
        s = 0.07
    elif n < 4000000:
        s = 0.05
    elif n < 5000000:
        s = 0.03
    else:
        s = 0.02
    return s


def sync_xylim_width(ax):
    """Sync x and y axis limit to make them have the same width."""
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xdelta = xmax - xmin
    ydelta = ymax - ymin
    flank = abs(xdelta - ydelta) / 2
    if xdelta > ydelta:
        ymin -= flank
        ymax += flank
        ax.set_ylim(ymin, ymax)
    else:
        xmin -= flank
        xmax += flank
        ax.set_xlim(xmin, xmax)
    return
