import plotly.express as px


def interactive_scatter(data, hue=None, coord_base='umap', continous_cmap='viridis', size=5):
    fig = px.scatter(data_frame=data,
                     x=f'{coord_base}_0',
                     y=f'{coord_base}_1',
                     color_continuous_scale=continous_cmap,
                     color=hue)
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      xaxis_showticklabels=False,
                      yaxis_showticklabels=False,
                      width=800,
                      height=800)
    fig.update_traces(marker_size=size)
    return fig
