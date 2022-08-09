:py:mod:`ALLCools.plot.genome_track.GtfTrack`
=============================================

.. py:module:: ALLCools.plot.genome_track.GtfTrack

.. autoapi-nested-parse::

   This GtfTrack class replaced the original pygenometracks class
   It supports read-in a gffutils database, rather than parse very time



Module Contents
---------------

.. py:data:: DEFAULT_BED_COLOR
   :annotation: = #1f78b4

   

.. py:data:: DISPLAY_BED_VALID
   :annotation: = ['collapsed', 'triangles', 'interleaved', 'stacked']

   

.. py:data:: DISPLAY_BED_SYNONYMOUS
   

   

.. py:data:: DEFAULT_DISPLAY_BED
   :annotation: = stacked

   

.. py:data:: AROUND_REGION
   :annotation: = 100000

   

.. py:function:: _is_sqlite3(path)


.. py:class:: ReadGtf(file_path, prefered_name='transcript_name', merge_transcripts=True)

   Bases: :py:obj:`object`

   Reads a gtf file.

   Example:
   gtf = ReadGtf("file.gtf")
   for interval in gtf:
       print interval.start


   .. py:method:: __iter__(self)


   .. py:method:: __next__(self)

      :return: bedInterval object


   .. py:method:: get_bed_interval(self)

      Process a transcript from the database,
      retrieve all the values and return
      a namedtuple object



.. py:class:: GtfTrack(*args, **kwarg)

   Bases: :py:obj:`pygenometracks.tracks.BedTrack.BedTrack`

   The GenomeTrack object is a holder for all tracks that are to be plotted.
   For example, to plot a bedgraph file a new class that extends GenomeTrack
   should be created.

   It is expected that all GenomeTrack objects have a plot method.


   .. py:attribute:: SUPPORTED_ENDINGS
      :annotation: = ['gtf', 'gtf.gz', 'gtf.db']

      

   .. py:attribute:: TRACK_TYPE
      :annotation: = gtf

      

   .. py:attribute:: OPTIONS_TXT
      

      

   .. py:attribute:: DEFAULTS_PROPERTIES
      

      

   .. py:attribute:: NECESSARY_PROPERTIES
      :annotation: = ['file']

      

   .. py:attribute:: SYNONYMOUS_PROPERTIES
      

      

   .. py:attribute:: POSSIBLE_PROPERTIES
      

      

   .. py:attribute:: BOOLEAN_PROPERTIES
      :annotation: = ['labels', 'merge_transcripts', 'global_max_row', 'arrowhead_included', 'all_labels_inside',...

      

   .. py:attribute:: STRING_PROPERTIES
      :annotation: = ['prefered_name', 'file', 'file_type', 'overlay_previous', 'orientation', 'title', 'style',...

      

   .. py:attribute:: FLOAT_PROPERTIES
      

      

   .. py:attribute:: INTEGER_PROPERTIES
      

      

   .. py:method:: set_properties_defaults(self)


   .. py:method:: get_bed_handler(self, plot_regions=None)



