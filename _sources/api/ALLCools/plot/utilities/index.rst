:py:mod:`ALLCools.plot.utilities`
=================================

.. py:module:: ALLCools.plot.utilities


Module Contents
---------------

.. py:function:: _density_based_sample(data: pandas.DataFrame, coords: list, portion=None, size=None, seed=None)

   Down sample data based on density, to prevent overplot in dense region and decrease plotting time.


.. py:function:: _translate_coord_name(coord_name)


.. py:function:: _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5)


.. py:function:: _extract_coords(data, coord_base, x, y)


.. py:function:: zoom_min_max(vmin, vmax, scale)

   Zoom min and max value.


.. py:function:: zoom_ax(ax, zoom_scale, on='both')

   Zoom ax on both x and y-axis.


.. py:function:: smart_number_format(x, pos=None)

   Format numbers automatically.


.. py:function:: add_ax_box(ax, expend=0, **patch_kws)

   Add a box around the ax.


.. py:function:: tight_hue_range(hue_data, portion)

   Automatic select a SMALLEST data range that covers [portion] of the data.


.. py:function:: _take_data_series(data, k)


.. py:function:: _auto_size(ax, n_dots)

   Auto determine dot size based on ax size and n dots


