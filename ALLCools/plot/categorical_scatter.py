import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from .color import level_one_palette
from .contour import density_contour
from .text_anno_scatter import _text_anno_scatter
from .utilities import (
    _auto_size,
    _density_based_sample,
    _extract_coords,
    _make_tiny_axis_label,
    _take_data_series,
    tight_hue_range,
    zoom_ax,
)


def categorical_scatter(
    data,
    ax=None,
    # coords
    coord_base="umap",
    x=None,
    y=None,
    # color
    hue=None,
    palette="auto",
    color=None,
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
    s="auto",
    size=None,
    sizes: dict = None,
    size_norm=None,
    size_portion=0.95,
    # other
    axis_format="tiny",
    max_points=50000,
    labelsize=4,
    linewidth=0,
    zoomxy=1.05,
    outline=None,
    outline_pad=3,
    outline_kws=None,
    scatter_kws=None,
    return_fig=False,
    rasterized="auto",
):
    """
    Plot categorical scatter plot with versatile options.

    Parameters
    ----------
    rasterized
        Whether to rasterize the figure.
    return_fig
        Whether to return the figure.
    size_portion
        The portion of the figure to be used for the size norm.
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
    hue : str
        categorical col name or series for color hue.
    palette : str or dict
        palette for color hue.
    color
        specify single color for all the dots
    text_anno
        categorical col name or series for text annotation.
    text_anno_kws
        kwargs for text annotation.
    text_anno_palette
        palette for text annotation.
    text_transform
        transform for text annotation.
    dodge_text
        whether to dodge text annotation.
    dodge_kws
        kwargs for dodge text annotation.
    show_legend
        whether to show legend.
    legend_kws
        kwargs for legend.
    s
        single size value of all the dots.
    size
        mappable size of the dots.
    sizes
        mapping size to the sizes value.
    size_norm
        normalize size range for mapping.
    axis_format
        axis format.
    max_points
        maximum number of points to plot.
    labelsize
        label size pass to `ax.text`
    linewidth
        line width pass to `ax.scatter`
    zoomxy
        zoom factor for x and y-axis.
    outline
        categorical col name or series for outline.
    outline_pad
        outline padding.
    outline_kws
        kwargs for outline.
    scatter_kws
        kwargs for scatter.

    Returns
    -------
    if return_fig is True, return the figure and axes.
    else, return None.
    """
    # init figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    else:
        fig = None

    # add coords
    _data, x, y = _extract_coords(data, coord_base, x, y)
    # _data has 2 cols: "x" and "y"

    # down sample plot data if needed.
    if max_points is not None:
        if _data.shape[0] > max_points:
            _data = _density_based_sample(_data, seed=1, size=max_points, coords=["x", "y"])
    n_dots = _data.shape[0]

    # determine rasterized
    if rasterized == "auto":
        if n_dots > 200:
            rasterized = True
        else:
            rasterized = False

    # auto size if user didn't provide one
    if s == "auto":
        s = _auto_size(ax, n_dots)

    # default scatter options
    _scatter_kws = {"linewidth": 0, "s": s, "legend": None, "palette": palette, "rasterized": rasterized}
    if color is not None:
        if hue is not None:
            raise ValueError("Only one of color and hue can be provided")
        _scatter_kws["color"] = color
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    # deal with color
    palette_dict = None
    if hue is not None:
        if isinstance(hue, str):
            _data["hue"] = _take_data_series(data, hue)
        else:
            _data["hue"] = hue.copy()
        hue = "hue"
        _data["hue"] = _data["hue"].astype("category").cat.remove_unused_categories()

        # if the object has get_palette method, use it (AnnotZarr)
        palette = _scatter_kws["palette"]
        if hasattr(_data, "get_palette"):
            palette_dict = _data.get_palette(hue)

        # deal with other color palette
        if palette_dict is None:
            if isinstance(palette, str) or isinstance(palette, list):
                palette_dict = level_one_palette(_data["hue"], order=None, palette=palette)
            elif isinstance(palette, dict):
                palette_dict = palette
            else:
                raise TypeError(f"Palette can only be str, list or dict, " f"got {type(palette)}")
        _scatter_kws["palette"] = palette_dict

    # deal with size
    if size is not None:
        if isinstance(size, str):
            _data["size"] = _take_data_series(data, size).astype(float)
        else:
            _data["size"] = size.astype(float)
        size = "size"

        if size_norm is None:
            # get the smallest range that include "size_portion" of data
            size_norm = tight_hue_range(_data["size"], size_portion)

            # snorm is the normalizer for size
            size_norm = Normalize(vmin=size_norm[0], vmax=size_norm[1])

        # discard s from _scatter_kws and use size in sns.scatterplot
        s = _scatter_kws.pop("s")
        if sizes is None:
            sizes = (min(s, 1), s)

    sns.scatterplot(
        x="x",
        y="y",
        data=_data,
        ax=ax,
        hue=hue,
        size=size,
        sizes=sizes,
        size_norm=size_norm,
        **_scatter_kws,
    )

    # deal with text annotation
    if text_anno is not None:
        if isinstance(text_anno, str):
            _data["text_anno"] = _take_data_series(data, text_anno)
        else:
            _data["text_anno"] = text_anno.copy()
        if str(_data["text_anno"].dtype) == "category":
            _data["text_anno"] = _data["text_anno"].cat.remove_unused_categories()

        _text_anno_scatter(
            data=_data[["x", "y", "text_anno"]],
            ax=ax,
            x="x",
            y="y",
            dodge_text=dodge_text,
            dodge_kws=dodge_kws,
            palette=text_anno_palette,
            text_transform=text_transform,
            anno_col="text_anno",
            text_anno_kws=text_anno_kws,
            labelsize=labelsize,
        )

    # deal with outline
    if outline:
        if isinstance(outline, str):
            _data["outline"] = _take_data_series(data, outline)
        else:
            _data["outline"] = outline.copy()
        _outline_kws = {
            "linewidth": linewidth,
            "palette": None,
            "c": "lightgray",
            "single_contour_pad": outline_pad,
        }
        if outline_kws is not None:
            _outline_kws.update(outline_kws)
        density_contour(ax=ax, data=_data, x="x", y="y", groupby="outline", **_outline_kws)

    # clean axis
    if axis_format == "tiny":
        _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=labelsize)
    elif (axis_format == "empty") or (axis_format is None):
        sns.despine(ax=ax, left=True, bottom=True)
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    else:
        pass

    # deal with legend
    if show_legend and (hue is not None):
        n_hue = len(palette_dict)
        _legend_kws = {
            "ncol": (1 if n_hue <= 20 else 2 if n_hue <= 40 else 3),
            "fontsize": labelsize,
            "bbox_to_anchor": (1.05, 1),
            "loc": "upper left",
            "borderaxespad": 0.0,
        }
        if legend_kws is not None:
            _legend_kws.update(legend_kws)

        handles = []
        labels = []
        exist_hues = _data["hue"].unique()
        for hue_name, color in palette_dict.items():
            if hue_name not in exist_hues:
                # skip hue_name that do not appear in the plot
                continue
            handle = Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=_legend_kws["fontsize"],
            )
            handles.append(handle)
            labels.append(hue_name)
        _legend_kws["handles"] = handles
        _legend_kws["labels"] = labels
        ax.legend(**_legend_kws)

    if zoomxy is not None:
        zoom_ax(ax, zoomxy)

    if return_fig:
        return (fig, ax), _data
    else:
        return
