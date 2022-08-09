:py:mod:`ALLCools.plot.genome_track.HiCMatrixCoolTrack`
=======================================================

.. py:module:: ALLCools.plot.genome_track.HiCMatrixCoolTrack


Module Contents
---------------

.. py:data:: DEFAULT_MATRIX_COLORMAP
   :annotation: = RdYlBu_r

   

.. py:data:: log
   

   

.. py:class:: HiCMatrixCoolTrack(*args, **kwargs)

   Bases: :py:obj:`pygenometracks.tracks.GenomeTrack.GenomeTrack`

   .. py:attribute:: SUPPORTED_ENDINGS
      :annotation: = ['.cool', '.mcool']

      

   .. py:attribute:: TRACK_TYPE
      :annotation: = cooler

      

   .. py:attribute:: OPTIONS_TXT
      

      

   .. py:attribute:: DEFAULTS_PROPERTIES
      

      

   .. py:attribute:: NECESSARY_PROPERTIES
      :annotation: = ['file']

      

   .. py:attribute:: SYNONYMOUS_PROPERTIES
      

      

   .. py:attribute:: POSSIBLE_PROPERTIES
      

      

   .. py:attribute:: BOOLEAN_PROPERTIES
      :annotation: = ['show_masked_bins', 'rasterize']

      

   .. py:attribute:: STRING_PROPERTIES
      :annotation: = ['file', 'file_type', 'overlay_previous', 'orientation', 'transform', 'title', 'colormap']

      

   .. py:attribute:: FLOAT_PROPERTIES
      

      

   .. py:attribute:: INTEGER_PROPERTIES
      

      

   .. py:method:: set_properties_defaults()


   .. py:method:: plot(ax, chrom_region, region_start, region_end)


   .. py:method:: plot_y_axis(cbar_ax, plot_ax)


   .. py:method:: pcolormesh_45deg(ax, matrix_c, start_pos_vector)

      Plot a 45 degree heatmap.

      Turns the matrix 45 degrees and adjusts the bins to match the actual start end positions.



