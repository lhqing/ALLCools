import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns
import copy
from matplotlib.cm import get_cmap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap, Normalize


def _continuous_color_palette(color, n, skip_border=1):
    """
    This function concatenate the result of both sns.light_palette
    and sns.dark_palette to get a wider color range
    """
    if n == 1:
        return [color]
    if n < 1:
        raise ValueError('parameter n colors must >= 1.')

    # this is just trying to make sure len(color) == n
    light_n = (n + 2 * skip_border) // 2
    light_colors = sns.light_palette(color, n_colors=light_n)[skip_border:]
    dark_n = n + 2 * skip_border - light_n + 1
    dark_colors = sns.dark_palette(color, n_colors=dark_n, reverse=True)[1:-skip_border]
    colors = light_colors + dark_colors
    return colors


def get_kv_dict(data_df, major, sub):
    _dict = data_df.groupby(major)[sub].unique().to_dict()
    return _dict


def level_one_palette(name_list, order=None, palette='auto'):
    name_set = set(name_list)
    if palette == 'auto':
        if len(name_set) < 10:
            palette = 'tab10'
        elif len(name_set) < 20:
            palette = 'tab20'
        else:
            palette = 'rainbow'

    if order is None:
        order = list(sorted(name_set))
    else:
        if (set(order) != name_set) or (len(order) != len(name_set)):
            raise ValueError('Order is not equal to set(name_list).')
    n = len(order)
    colors = sns.color_palette(palette, n)
    color_palette = {}
    for name, color in zip(order, colors):
        color_palette[name] = color
    return color_palette


def level_two_palette(major_color, major_sub_dict,
                      major_order=None, palette='auto',
                      skip_border_color=2):
    if isinstance(major_color, list):
        major_color_dict = level_one_palette(major_color, palette=palette, order=major_order)
    else:
        major_color_dict = major_color

    sub_id_list = []
    for subs in major_sub_dict.values():
        sub_id_list += list(subs)
    if len(sub_id_list) != len(set(sub_id_list)):
        raise ValueError('Sub id in the major_dub_dict is not unique.')

    color_palette = {}
    for major, color in major_color_dict.items():
        subs = major_sub_dict[major]
        n = len(subs)
        colors = _continuous_color_palette(color, n, skip_border=skip_border_color)
        for sub, _color in zip(subs, colors):
            color_palette[sub] = _color
    return color_palette


def palplot(pal, transpose=False):
    if transpose:
        fig, ax = plt.subplots(figsize=(1, len(pal)))
    else:
        fig, ax = plt.subplots(figsize=(len(pal), 1))

    plot_color_legend(pal, ax, order=None, transpose=transpose)
    return fig, ax


def plot_colorbar(cax, cmap, hue_norm, cnorm=None, label=None, orientation='vertical',
                  labelsize=4, linewidth=0.5):
    if isinstance(cmap, str):
        cmap = copy.copy(get_cmap(cmap))
    if cnorm is None:
        cnorm = Normalize(vmin=hue_norm[0],
                          vmax=hue_norm[1])
    from .utilities import smart_number_format

    colorbar = ColorbarBase(cax,
                            cmap=cmap,
                            norm=cnorm,
                            format=ticker.FuncFormatter(smart_number_format),
                            orientation=orientation,
                            extend='both')
    colorbar.locator = ticker.MaxNLocator(nbins=3)
    colorbar.update_ticks()

    colorbar.set_label(label, fontsize=labelsize)
    colorbar.outline.set_linewidth(linewidth)
    colorbar.ax.tick_params(size=labelsize,
                            labelsize=labelsize,
                            width=linewidth)
    return cax


def plot_color_legend(palette, ax, order=None, interpolation=None, transpose=False):
    if order is None:
        order = palette.keys()

    colors = []
    for c in order:
        colors.append(palette[c])

    n = len(order)
    data = np.arange(n).reshape(1, n)
    if transpose:
        data = data.T
    ax.imshow(data, interpolation=interpolation, aspect="auto",
              cmap=ListedColormap(colors))
    if not transpose:
        ax.set(xticklabels=list(order),
               xticks=range(0, n),
               yticks=[])
        ax.xaxis.set_tick_params(labelrotation=90)
    else:
        ax.set(yticklabels=list(order),
               yticks=range(0, n),
               xticks=[])
    return
