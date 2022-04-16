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

   The GenomeTrack object is a holder for all tracks that are to be plotted.
   For example, to plot a bedgraph file a new class that extends GenomeTrack
   should be created.

   It is expected that all GenomeTrack objects have a plot method.


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
      

      

   .. py:method:: set_properties_defaults(self)


   .. py:method:: plot(self, ax, chrom_region, region_start, region_end)


   .. py:method:: plot_y_axis(self, cbar_ax, plot_ax)

      Plot the scale of the y axis with respect to the plot_axis
      :param ax: axis to use to plot the scale
      :param plot_axis: the reference axis to get the max and min.
      :param transform: what was the transformation of the data
      :param log_pseudocount:
      :param y_axis: 'tranformed' or 'original'
      :param only_at_ticks: False: only min_max are diplayed
                            True: only ticks values are displayed

      Returns:



   .. py:method:: pcolormesh_45deg(self, ax, matrix_c, start_pos_vector)

      Turns the matrix 45 degrees and adjusts the
      bins to match the actual start end positions.



