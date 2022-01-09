:py:mod:`ALLCools.igv.Track`
============================

.. py:module:: ALLCools.igv.Track


Module Contents
---------------

.. py:class:: DataTrack(**kwargs)

   Bases: :py:obj:`ALLCools.igv.BaseClass.Track`

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


   .. py:method:: add_data_range(self, **kwargs)



.. py:class:: FeatureTrack(**kwargs)

   Bases: :py:obj:`ALLCools.igv.BaseClass.Track`

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



