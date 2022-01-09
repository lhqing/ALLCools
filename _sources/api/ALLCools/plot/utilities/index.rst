:py:mod:`ALLCools.plot.utilities`
=================================

.. py:module:: ALLCools.plot.utilities


Module Contents
---------------

.. py:function:: _density_based_sample(data: pandas.DataFrame, coords: list, portion=None, size=None, seed=None)

   down sample data based on density, to prevent overplot in dense region and decrease plotting time


.. py:function:: _translate_coord_name(coord_name)


.. py:function:: _make_tiny_axis_label(ax, x, y, arrow_kws=None, fontsize=5)

   this function assume coord is [0, 1]


.. py:function:: _extract_coords(data, coord_base, x, y)


.. py:function:: zoom_min_max(vmin, vmax, scale)


.. py:function:: zoom_ax(ax, zoom_scale, on='both')


.. py:function:: smart_number_format(x, pos=None)


.. py:function:: add_ax_box(ax, expend=0, **patch_kws)


.. py:function:: tight_hue_range(hue_data, portion)

   Automatic select a SMALLEST data range that covers [portion] of the data


