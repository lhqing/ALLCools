:py:mod:`ALLCools.plot.qc_plots`
================================

.. py:module:: ALLCools.plot.qc_plots


Module Contents
---------------

.. py:function:: plot_on_plate(data, value_col, groupby, ncols=4, plate_base=384, figsize=(5, 3.5), row_base='Row384', col_base='Col384', vmin=0, vmax=1, heatmap_kws=None, aggregation_func=None)

   Plot metadata into 384 or 96 plate view (heatmap)
   :param data: dataframe contain all kinds of metadata
   :param value_col: value to be plotted on plate view
   :param groupby: groupby column, typically groupby plate id column(s) to plot each plate separately
   :param ncols: number of column for axes, nrows will be calculated accordingly
   :param plate_base: {384, 96} size of the plate view
   :param figsize: matplotlib.Figure figsize
   :param vmin: cmap vmin
   :param vmax: cmap vmax
   :param heatmap_kws: kws pass to sns.heatmap
   :param aggregation_func: apply to reduce rows after groupby if the row is not unique


.. py:function:: simple_violin(data, x, y, rotate_x=True)


.. py:function:: cutoff_vs_cell_remain(data, name='', xlim_quantile=(0.001, 0.999), ylim=None, bins=100)


.. py:function:: success_vs_fail(data, filter_col, filter_cutoff, x, y, ax)


.. py:function:: plot_dispersion(data, hue='gene_subset', zlab='dispersion', data_quantile=(0.01, 0.99), save_animate_path=None, fig_kws=None)


.. py:function:: plot_hvf_selection()


