import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


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


def level_one_palette(name_list, order=None, palette='default'):
    name_set = set(name_list)
    if palette == 'default':
        if len(set(name_list)) < 10:
            palette = 'tab10'
        else:
            palette = 'tab20'

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
                      major_order=None, palette='default',
                      skip_border_color=2):
    if isinstance(major_color, list):
        if len(major_color) > 20:
            print(f'Warning: too much major color {len(major_color)} make the palette less distinguishable.')
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
    n = len(pal)
    data = np.arange(n).reshape(1, n)
    if transpose:
        data = data.T
    ax.imshow(data, interpolation="nearest", aspect="auto",
              cmap=mpl.colors.ListedColormap(list(pal.values())))
    if not transpose:
        ax.set(xticklabels=list(pal.keys()),
               xticks=range(0, len(pal)),
               yticks=[])
        ax.xaxis.set_tick_params(labelrotation=90)
    else:
        ax.set(yticklabels=list(pal.keys()),
               yticks=range(0, len(pal)),
               xticks=[])
    return fig, ax
