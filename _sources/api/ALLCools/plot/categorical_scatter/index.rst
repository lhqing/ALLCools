:py:mod:`ALLCools.plot.categorical_scatter`
===========================================

.. py:module:: ALLCools.plot.categorical_scatter


Module Contents
---------------

.. py:function:: categorical_scatter(data, ax, coord_base='umap', x=None, y=None, hue=None, palette='auto', text_anno=None, text_anno_kws=None, text_anno_palette=None, text_transform=None, dodge_text=False, dodge_kws=None, show_legend=False, legend_kws=None, s=5, size=None, sizes: dict = None, size_norm=None, axis_format='tiny', max_points=5000, labelsize=4, linewidth=0, zoomxy=1.05, outline=None, outline_pad=3, outline_kws=None, scatter_kws=None, return_fig=False)

   Plot scatter plot with these options:
       - Color by a categorical variable, and generate legend of the variable if needed
       - Add text annotation using a categorical variable
       - Circle categories with outlines

   :param data: Dataframe that contains coordinates and categorical variables
   :param ax: this function do not generate ax, must provide an ax
   :param coord_base: coords name, if provided, will automatically search for x and y
   :param x: x coord name
   :param y: y coord name
   :param hue: categorical col name or series for color hue
   :param palette: palette of the hue, str or dict
   :param text_anno: categorical col name or series for text annotation
   :param text_anno_kws:
   :param text_anno_palette:
   :param text_transform:
   :param dodge_text:
   :param dodge_kws:
   :param show_legend:
   :param legend_kws:
   :param s:
   :param size:
   :param sizes:
   :param size_norm:
   :param axis_format:
   :param max_points:
   :param labelsize:
   :param linewidth:
   :param zoomxy:
   :param outline:
   :param outline_pad:
   :param outline_kws:
   :param scatter_kws: kws dict pass to sns.scatterplot


