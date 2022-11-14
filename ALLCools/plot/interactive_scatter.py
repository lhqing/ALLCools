import anndata


def interactive_scatter(
    data, hue=None, coord_base="tsne", continuous_cmap="viridis", size=5, max_points=3000, hover_data=("obs_names",)
):
    """
    Plot an interactive scatter plot with plotly.

    Parameters
    ----------
    data :
        AnnData object or pd.DataFrame
    hue :
        Column name in data to color the points
    coord_base :
        Coordinate base to plot. Default is tsne.
    continuous_cmap :
        Continuous colormap to use. Default is viridis.
    size :
        Size of the points. Default is 5.
    max_points :
        Maximum number of points to plot. Default is 3000.
    hover_data :
        Column names in data to show in the hover data. Default is (obs_names, ).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    try:
        import plotly.express as px
    except ImportError:
        raise ImportError("Please install plotly to use this function.")
    if isinstance(data, anndata.AnnData):
        _data = data.obs.copy()
        for k in hover_data:
            if k == "obs_names":
                _data["obs_names"] = data.obs_names
            else:
                _data[k] = data.obs[k]
        _data[f"{coord_base}_0"] = data.obsm[f"X_{coord_base}"][:, 0]
        _data[f"{coord_base}_1"] = data.obsm[f"X_{coord_base}"][:, 1]
    else:
        _data = data.copy()
    if isinstance(max_points, int) and (max_points < _data.shape[0]):
        _data = _data.sample(max_points)

    fig = px.scatter(
        data_frame=_data,
        x=f"{coord_base}_0",
        y=f"{coord_base}_1",
        color_continuous_scale=continuous_cmap,
        color=hue,
        hover_data=hover_data,
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
