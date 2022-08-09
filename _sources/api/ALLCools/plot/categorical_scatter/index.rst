:py:mod:`ALLCools.plot.categorical_scatter`
===========================================

.. py:module:: ALLCools.plot.categorical_scatter


Module Contents
---------------

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


