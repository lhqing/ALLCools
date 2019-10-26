import pandas as pd


def _text_anno_scatter(data: pd.DataFrame,
                       ax,
                       x: str,
                       y: str,
                       edge_color=(0.5, 0.5, 0.5, 0.2),
                       face_color=(0.8, 0.8, 0.8, 0.2),
                       palette: dict = None,
                       dodge_text=False,
                       anno_col='text_anno',
                       text_anno_kws=None,
                       text_transform=None,
                       dodge_kws=None,
                       linewidth=0.5,
                       labelsize=5):
    """Add text annotation to a scatter plot"""
    # prepare kws
    _text_anno_kws = dict(fontsize=labelsize,
                          fontweight='black',
                          horizontalalignment='center',
                          verticalalignment='center')
    if text_anno_kws is not None:
        _text_anno_kws.update(text_anno_kws)

    # plot each text
    text_list = []
    for text, sub_df in data.groupby(anno_col):
        if text_transform is None:
            text = str(text)
        else:
            text = text_transform(text)
        if text.lower() in ['', 'nan']:
            continue
        _x, _y = sub_df[[x, y]].median()

        if palette is not None:
            _fc = palette[text]
        else:
            _fc = face_color
        text = ax.text(_x, _y, text,
                       fontdict=_text_anno_kws,
                       bbox=dict(boxstyle="round",
                                 ec=edge_color,
                                 fc=_fc,
                                 linewidth=linewidth))
        text_list.append(text)

    if dodge_text:
        try:
            from adjustText import adjust_text

            _dodge_parms = dict(force_points=(0.02, 0.05),
                                arrowprops=dict(arrowstyle="->",
                                                fc=edge_color,
                                                ec="none",
                                                connectionstyle="angle,angleA=-90,angleB=180,rad=5"),
                                autoalign='xy')
            if dodge_kws is not None:
                _dodge_parms.update(dodge_kws)
            adjust_text(text_list, x=data['x'], y=data['y'], **_dodge_parms)
        except ModuleNotFoundError:
            print('Install adjustText package to dodge text, see its github page for help')

    return
