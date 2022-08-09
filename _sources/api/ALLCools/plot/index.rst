:py:mod:`ALLCools.plot`
=======================

.. py:module:: ALLCools.plot


Subpackages
-----------
.. toctree::
   :titlesonly:
   :maxdepth: 3

   genome_track/index.rst


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   categorical_scatter/index.rst
   color/index.rst
   continuous_scatter/index.rst
   contour/index.rst
   decomposition/index.rst
   dendro/index.rst
   interactive_scatter/index.rst
   qc_plots/index.rst
   size/index.rst
   sunburst/index.rst
   text_anno_scatter/index.rst
   utilities/index.rst


Package Contents
----------------

.. py:function:: categorical_scatter(data, ax=None, coord_base='umap', x=None, y=None, hue=None, palette='auto', color=None, text_anno=None, text_anno_kws=None, text_anno_palette=None, text_transform=None, dodge_text=False, dodge_kws=None, show_legend=False, legend_kws=None, s='auto', size=None, sizes: dict = None, size_norm=None, size_portion=0.95, axis_format='tiny', max_points=50000, labelsize=4, linewidth=0, zoomxy=1.05, outline=None, outline_pad=3, outline_kws=None, scatter_kws=None, return_fig=False, rasterized='auto')

   Plot categorical scatter plot with versatile options.

   :param rasterized: Whether to rasterize the figure.
   :param return_fig: Whether to return the figure.
   :param size_portion: The portion of the figure to be used for the size norm.
   :param data: Dataframe that contains coordinates and categorical variables
   :param ax: this function do not generate ax, must provide an ax
   :param coord_base: coords name, if provided, will automatically search for x and y
   :param x: x coord name
   :param y: y coord name
   :param hue: categorical col name or series for color hue.
   :type hue: str
   :param palette: palette for color hue.
   :type palette: str or dict
   :param color: specify single color for all the dots
   :param text_anno: categorical col name or series for text annotation.
   :param text_anno_kws: kwargs for text annotation.
   :param text_anno_palette: palette for text annotation.
   :param text_transform: transform for text annotation.
   :param dodge_text: whether to dodge text annotation.
   :param dodge_kws: kwargs for dodge text annotation.
   :param show_legend: whether to show legend.
   :param legend_kws: kwargs for legend.
   :param s: single size value of all the dots.
   :param size: mappable size of the dots.
   :param sizes: mapping size to the sizes value.
   :param size_norm: normalize size range for mapping.
   :param axis_format: axis format.
   :param max_points: maximum number of points to plot.
   :param labelsize: label size pass to `ax.text`
   :param linewidth: line width pass to `ax.scatter`
   :param zoomxy: zoom factor for x and y-axis.
   :param outline: categorical col name or series for outline.
   :param outline_pad: outline padding.
   :param outline_kws: kwargs for outline.
   :param scatter_kws: kwargs for scatter.

   :returns: * *if return_fig is True, return the figure and axes.*
             * *else, return None.*


.. py:function:: continuous_scatter(data, ax=None, coord_base='umap', x=None, y=None, scatter_kws=None, hue=None, hue_norm=None, hue_portion=0.95, color=None, cmap='viridis', colorbar=True, colorbar_label_kws=None, size=None, size_norm=None, size_portion=0.95, sizes=None, sizebar=True, text_anno=None, dodge_text=False, dodge_kws=None, text_anno_kws=None, text_anno_palette=None, text_transform=None, axis_format='tiny', max_points=50000, s='auto', labelsize=4, linewidth=0.5, cax=None, zoomxy=1.05, outline=None, outline_kws=None, outline_pad=2, return_fig=False, rasterized='auto')


.. py:function:: plot_decomp_scatters(adata, n_components, base_name='PC', obsm_name='X_pca', hue=None, palette='viridis', hue_quantile=(0.25, 0.75), nrows=5, ncols=5)


.. py:function:: plot_dendrogram(linkage_df, ax, dendro=None, labels=None, dendro_kws=None, plot_node_id=False, plot_non_singleton=True, plot_kws=None, node_hue=None, node_hue_norm=None, node_hue_cbar=True, node_hue_cbar_frac=0.1, node_palette='viridis', node_size=None, node_size_norm=None, node_sizes=None, line_hue=None, line_hue_norm=None, line_palette='gray_r', linewidth=1.5, edge_color='gray', marker_size=60, marker_color='lightblue')

   Plot a dendrogram with customizations.

   :param linkage_df:
   :param dendro:
   :param labels:
   :param dendro_kws:
   :param ax:
   :param plot_node_id:
   :param plot_non_singleton:
   :param plot_kws:
   :param node_hue:
   :param node_hue_norm:
   :param node_hue_cbar:
   :param node_hue_cbar_frac:
   :param node_palette:
   :param node_size:
   :param node_size_norm:
   :param node_sizes:
   :param line_hue:
   :param line_hue_norm:
   :param line_palette:
   :param linewidth:
   :param edge_color:
   :param marker_size:
   :param marker_color:

   :returns: * *dendro* -- dendrogram object
             * *node_pos* -- node position


.. py:function:: interactive_scatter(data, hue=None, coord_base='umap', continous_cmap='viridis', size=5, max_points=3000)

   Plot an interactive scatter plot with plotly.


.. py:function:: sunbrust(pie_data, ax, hue=None, hue_portion=0.5, cmap='viridis', colorbar=True, colorbar_kws=None, inner_radius=0.25, outer_radius=1, anno_col=None, text_anno='text', anno_layer_size=0.05, col_color_dict=None, startangle=0, anno_ang_min=5, anno_border=1.2, text_expend=1.05, uniform_section=False, order_dict=None)

   Plot a sunburst diagram.


.. py:function:: add_ax_box(ax, expend=0, **patch_kws)

   Add a box around the ax.


