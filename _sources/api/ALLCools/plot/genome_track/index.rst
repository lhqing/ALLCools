:py:mod:`ALLCools.plot.genome_track`
====================================

.. py:module:: ALLCools.plot.genome_track


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   _pygenometrack/index.rst
   utilities/index.rst


Package Contents
----------------

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



.. py:function:: read_gtf(gtf_path)


.. py:function:: subset_gtf(gtf, regions, output_path=None, select_feature=None)


