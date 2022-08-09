:py:mod:`ALLCools.plot.genome_track._pygenometrack`
===================================================

.. py:module:: ALLCools.plot.genome_track._pygenometrack

.. autoapi-nested-parse::

   Modified from pygenometracks
   LICENSE: https://github.com/deeptools/pyGenomeTracks/blob/master/LICENSE



Module Contents
---------------

.. py:data:: FORMAT
   :annotation: = [%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s

   

.. py:data:: log
   

   

.. py:data:: DEFAULT_TRACK_HEIGHT
   :annotation: = 0.5

   

.. py:data:: DEFAULT_FIGURE_WIDTH
   :annotation: = 40

   

.. py:data:: DEFAULT_WIDTH_RATIOS
   :annotation: = [0.01, 0.9, 0.1]

   

.. py:data:: DEFAULT_MARGINS
   

   

.. py:class:: MultiDict

   Bases: :py:obj:`collections.OrderedDict`

   Class to allow identically named
   sections in configuration file
   by appending the section number
   for example:
   1. section name

   .. py:attribute:: _unique
      :annotation: = 0

      

   .. py:method:: __setitem__(self, key, val)

      Set self[key] to value.



.. py:class:: PlotTracks(config, fig_width=DEFAULT_FIGURE_WIDTH, fig_height=None, fontsize=None, dpi=None, track_label_width=None, plot_regions=None, plot_width=None)

   Bases: :py:obj:`object`

   .. py:method:: get_available_tracks()
      :staticmethod:


   .. py:method:: get_tracks_height(self, start_region=None, end_region=None)

      The main purpose of the following loop is
      to get the height of each of the tracks
      because for the Hi-C the height is variable with respect
      to the range being plotted, the function is called
      when each plot is going to be printed.

      :param start_region: start of the region to plot.
                           Only used in case the plot is a Hi-C matrix
      :param end_region: end of the region to plot.
                         Only used in case the plot is a Hi-C matrix

      Returns:



   .. py:method:: plot(self, chrom, start, end, file_name=None, title=None, h_align_titles='left', decreasing_x_axis=False)


   .. py:method:: plot_vlines(self, axis_list, chrom_region, start_region, end_region)

      Plots dotted lines from the top of the first plot to the bottom
      of the last plot at the specified positions.

      :param axis_list: list of plotted axis
      :param chrom_region chromosome name
      :param start_region start position
      :param end_region end position

      :return: None


   .. py:method:: parse_tracks(self, tracks_file, plot_regions=None)

      Parses a configuration file

      :param tracks_file: file path containing the track configuration
      :param plot_regions: a list of tuple [(chrom1, start1, end1), (chrom2, start2, end2)]
                           on which the data should be loaded
                           here the vlines
      :return: array of dictionaries and vlines_file.
               One dictionary per track


   .. py:method:: close_files(self)

      Close all opened files


   .. py:method:: check_file_exists(track_dict, tracks_path, is_hic=False)
      :staticmethod:

      Checks if a file or list of files exists. If the file does not exist
      tries to check if the file may be relative to the track_file path,
      in such case the path is updated.
      :param track_dict: dictionary of track values. Should contain
                          a 'file' key containing the path of the file
                          or files to be checked separated by space
                          For example: file1 file2 file3
      :param tracks_path: path of the tracks file
      :param is_hic:
      :return: None


   .. py:method:: guess_filetype(track_dict, available_tracks)
      :staticmethod:

      :param track_dict: dictionary of track values with the 'file' key
                  containing a string path of the file or files.
                  Only the ending of the last file is used
                  in case when there are more files
      :param available_tracks: list of available tracks

      :return: string file type detected


   .. py:method:: cm2inch(*tupl)
      :staticmethod:


   .. py:method:: print_elapsed(start)
      :staticmethod:



.. py:class:: SpacerTrack(properties_dict)

   Bases: :py:obj:`pygenometracks.tracks.GenomeTrack.GenomeTrack`

   The GenomeTrack object is a holder for all tracks that are to be plotted.
   For example, to plot a bedgraph file a new class that extends GenomeTrack
   should be created.

   It is expected that all GenomeTrack objects have a plot method.


   .. py:attribute:: SUPPORTED_ENDINGS
      :annotation: = []

      

   .. py:attribute:: TRACK_TYPE
      :annotation: = spacer

      

   .. py:attribute:: DEFAULTS_PROPERTIES
      

      

   .. py:attribute:: NECESSARY_PROPERTIES
      :annotation: = []

      

   .. py:attribute:: SYNONYMOUS_PROPERTIES
      

      

   .. py:attribute:: POSSIBLE_PROPERTIES
      

      

   .. py:attribute:: BOOLEAN_PROPERTIES
      :annotation: = []

      

   .. py:attribute:: STRING_PROPERTIES
      :annotation: = ['overlay_previous', 'title', 'file_type']

      

   .. py:attribute:: FLOAT_PROPERTIES
      

      

   .. py:attribute:: INTEGER_PROPERTIES
      

      

   .. py:method:: plot(self, ax, chrom_region, start_region, end_region)


   .. py:method:: plot_y_axis(self, ax, plot_ax)

      Plot the scale of the y axis with respect to the plot_axis
      :param ax: axis to use to plot the scale
      :param plot_axis: the reference axis to get the max and min.
      :param transform: what was the transformation of the data
      :param log_pseudocount:
      :param y_axis: 'tranformed' or 'original'
      :param only_at_ticks: False: only min_max are diplayed
                            True: only ticks values are displayed

      Returns:




.. py:class:: XAxisTrack(*args, **kwargs)

   Bases: :py:obj:`pygenometracks.tracks.GenomeTrack.GenomeTrack`

   The GenomeTrack object is a holder for all tracks that are to be plotted.
   For example, to plot a bedgraph file a new class that extends GenomeTrack
   should be created.

   It is expected that all GenomeTrack objects have a plot method.


   .. py:attribute:: SUPPORTED_ENDINGS
      :annotation: = []

      

   .. py:attribute:: TRACK_TYPE
      :annotation: = x_axis

      

   .. py:attribute:: NECESSARY_PROPERTIES
      :annotation: = []

      

   .. py:attribute:: DEFAULTS_PROPERTIES
      

      

   .. py:attribute:: SYNONYMOUS_PROPERTIES
      

      

   .. py:attribute:: POSSIBLE_PROPERTIES
      

      

   .. py:attribute:: BOOLEAN_PROPERTIES
      :annotation: = []

      

   .. py:attribute:: STRING_PROPERTIES
      :annotation: = ['overlay_previous', 'title', 'where', 'file_type']

      

   .. py:attribute:: FLOAT_PROPERTIES
      

      

   .. py:attribute:: INTEGER_PROPERTIES
      

      

   .. py:method:: plot(self, ax, chrom_region, region_start, region_end)


   .. py:method:: plot_y_axis(self, ax, plot_ax)

      Plot the scale of the y axis with respect to the plot_axis
      :param ax: axis to use to plot the scale
      :param plot_axis: the reference axis to get the max and min.
      :param transform: what was the transformation of the data
      :param log_pseudocount:
      :param y_axis: 'tranformed' or 'original'
      :param only_at_ticks: False: only min_max are diplayed
                            True: only ticks values are displayed

      Returns:




.. py:function:: prepare_config(file_configs, add_spacer=True, spacer_height=0.5)

   Prepare pyGenomeTracks config string.

   :param file_configs: A list of file paths, or a list of file config dicts, containing file path and other config.
                        See pyGenomeTracks documentation for possible parameters.
   :param add_spacer: Whether add spacer between tracks
   :param spacer_height: spacer height in cm

   :returns:
   :rtype: a single config string that can be read by PlotTracks


.. py:function:: _prepare_track_config(file_config, track_type, track_class, add_spacer=True, spacer_height=0.5)


