:py:mod:`ALLCools.sandbox.igv.Session`
======================================

.. py:module:: ALLCools.sandbox.igv.Session


Module Contents
---------------

.. py:function:: get_panel_int(panel)


.. py:class:: Session(genome='mm10', hasGeneTrack='false', hasSequenceTrack='true', locus='chr1:19000000-19100000', version='8', refseq_track=True)

   Bases: :py:obj:`ALLCools.sandbox.igv.BaseClass.Element`

   An XML element.

   This class is the reference implementation of the Element interface.

   An element's length is its number of subelements.  That means if you
   want to check if an element is truly empty, you should check BOTH
   its length AND its text attribute.

   The element tag, attribute names, and attribute values can be either
   bytes or strings.

   *tag* is the element name.  *attrib* is an optional dictionary containing
   element attributes. *extra* are additional element attributes given as
   keyword arguments.

   Example form:
       <tag attrib>text<child/>...</tag>tail


   .. py:method:: __add_resources(self)


   .. py:method:: __add_data_panel(self)


   .. py:method:: __add_feature_panel(self)


   .. py:method:: __add_other_routine(self)


   .. py:method:: __add_resource(self, path)


   .. py:method:: __add_sequence_track(self)


   .. py:method:: __add_refseq_gene(self)


   .. py:method:: add_data_track(self, track_kws, data_range_kws, panel='Data')


   .. py:method:: add_feature_track(self, track_kws, panel='Feature')



