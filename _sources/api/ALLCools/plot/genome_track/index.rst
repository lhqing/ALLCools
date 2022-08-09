:py:mod:`ALLCools.plot.genome_track`
====================================

.. py:module:: ALLCools.plot.genome_track


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   GtfTrack/index.rst
   HiCMatrixCoolTrack/index.rst
   _pygenometrack/index.rst


Package Contents
----------------

.. py:class:: PlotTracks(config, fig_width=DEFAULT_FIGURE_WIDTH, fig_height=None, fontsize=None, dpi=None, track_label_width=None, plot_regions=None, plot_width=None)

   Plotting functions.

   .. py:method:: get_available_tracks()
      :staticmethod:

      Get the available track types.


   .. py:method:: _get_tracks_height(start_region=None, end_region=None)

      Get track height.

      The main purpose of the following loop is to get the height of each of the tracks
      because for the Hi-C the height is variable with respect
      to the range being plotted, the function is called
      when each plot is going to be printed.

      :param start_region: Start of the region to plot. Only used in case the plot is a Hi-C matrix
      :param end_region: End of the region to plot. Only used in case the plot is a Hi-C matrix

      :returns: **track_height**
      :rtype: list


   .. py:method:: plot(chrom, start, end, file_name=None, title=None, h_align_titles='left', decreasing_x_axis=False)

      Make the genome track plot.


   .. py:method:: plot_vlines(axis_list, chrom_region, start_region, end_region)

      Plot vertical lines on the plot.

      Plots dotted lines from the top of the first plot to the bottom
      of the last plot at the specified positions.

      :param axis_list: list of plotted axis
      :param chrom_region: chromosome name
      :param start_region: start position
      :param end_region: end position


   .. py:method:: parse_tracks(tracks_file, plot_regions=None)

      Parse a configuration file.

      :param tracks_file: file path containing the track configuration
      :param plot_regions: a list of tuple [(chrom1, start1, end1), (chrom2, start2, end2)]
                           on which the data should be loaded here the vlines


   .. py:method:: close_files()

      Close all opened files.


   .. py:method:: check_file_exists(track_dict, tracks_path, is_hic=False)
      :staticmethod:

      Check if a file or list of files exists.

      If the file does not exist tries to check if the file may be relative to the track_file path,
      in such case the path is updated.

      :param track_dict: dictionary of track values. Should contain a 'file' key containing the path of the file
                         or files to be checked separated by space. For example: file1 file2 file3
      :param tracks_path: path of the tracks file
      :param is_hic: boolean indicating if the file is a hic matrix


   .. py:method:: guess_filetype(track_dict, available_tracks)
      :staticmethod:

      Guess the file type of the track.

      :param track_dict: dictionary of track values with the 'file' key containing a string path of the file or files.
                         Only the ending of the last file is used in case when there are more files
      :param available_tracks: list of available tracks

      :rtype: string file type detected


   .. py:method:: cm2inch(*tupl)
      :staticmethod:

      CM to INCH conversion.


   .. py:method:: print_elapsed(start)
      :staticmethod:

      Print elapsed time.



.. py:function:: prepare_config(file_configs, add_spacer=True, spacer_height=0.5)

   Prepare pyGenomeTracks config string.

   :param file_configs: A list of file paths, or a list of file config dicts, containing file path and other config.
                        See pyGenomeTracks documentation for possible parameters.
   :param add_spacer: Whether add spacer between tracks
   :param spacer_height: spacer height in cm

   :rtype: a single config string that can be read by PlotTracks


