:py:mod:`ALLCools.igv.BaseClass`
================================

.. py:module:: ALLCools.igv.BaseClass

.. autoapi-nested-parse::

   See IGV doc about xml, but this doc is actually out-of-date...
   https://software.broadinstitute.org/software/igv/Sessions

   I learn this by save IGV session xml file from IGV 2.17.2
   These session file is not in a stable API, use it cross IGV version do cause problems...



Module Contents
---------------

.. py:class:: Resources

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: Resource(path)

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: Panel(**kwargs)

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: Track(**kwargs)

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: DataRange(**kwargs)

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: PanelLayout(dividerFractions=0.9)

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



.. py:class:: HiddenAttributes

   Bases: :py:obj:`xml.etree.ElementTree.Element`

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



