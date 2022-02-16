import anndata
import plotly.express as px


def interactive_scatter(
        data,
        hue=None,
        coord_base="umap",
        continous_cmap="viridis",
        size=5,
        max_points=3000
):
    """
    Plot an interactive scatter plot with plotly

    Parameters
    ----------
    data
    hue
    coord_base
    continous_cmap
    size
    max_points

    Returns
    -------

    """
    if isinstance(data, anndata.AnnData):
        _data = data.obs.copy()
        _data[f'{coord_base}_0'] = data.obsm[f'X_{coord_base}'][:, 0]
        _data[f'{coord_base}_1'] = data.obsm[f'X_{coord_base}'][:, 1]
    else:
        _data = data.copy()
    if isinstance(max_points, int) and (max_points < _data.shape[0]):
        _data = _data.sample(max_points)

    fig = px.scatter(
        data_frame=_data,
        x=f'{coord_base}_0',
        y=f'{coord_base}_1',
        color_continuous_scale=continous_cmap,
        color=hue,
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_showticklabels=False,
        yaxis_showticklabels=False,
        width=800,
        height=800,
    )
    fig.update_traces(marker_size=size)
    return fig
