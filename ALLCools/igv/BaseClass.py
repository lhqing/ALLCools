"""
See IGV doc about xml, but this doc is actually out-of-date...
https://software.broadinstitute.org/software/igv/Sessions

I learn this by save IGV session xml file from IGV 2.17.2
These session file is not in a stable API, use it cross IGV version do cause problems...
"""

from xml.etree.ElementTree import Element
from .defaults import *


class Resources(Element):
    def __init__(self):
        _tag = 'Resources'
        super().__init__(_tag, attrib={})


class Resource(Element):
    def __init__(self, path):
        _tag = 'Resource'
        _attrib = {'path': path}
        super().__init__(_tag, attrib=_attrib)


class Panel(Element):
    def __init__(self,
                 **kwargs):
        _tag = 'Panel'
        super().__init__(_tag, attrib=kwargs)


class Track(Element):
    def __init__(self,
                 **kwargs):
        _tag = 'Track'
        for k in TRACK_REQUIRED_KEYS:
            if k not in kwargs:
                raise KeyError(f'Track missing argument {k}')

        # Default track parameters
        super().__init__(_tag, attrib=kwargs)


class DataRange(Element):
    def __init__(self,
                 **kwargs):
        _tag = 'DataRange'
        # Default DataRange parameters
        _attrib = DATARANGE_DEFAULT.copy()
        _attrib.update(kwargs)
        super().__init__(_tag, attrib=_attrib)


class PanelLayout(Element):
    def __init__(self,
                 dividerFractions=0.9):
        _tag = 'PanelLayout'
        # Default track parameters
        _attrib = dict(
            dividerFractions=str(dividerFractions))
        super().__init__(_tag, attrib=_attrib)


class HiddenAttributes(Element):
    def __init__(self):
        _tag = 'HiddenAttributes'
        super().__init__(_tag)
        # I don't know what IGV is doing here, but seems only stable lines so just copy here
        self.append(Element('Attribute', name="DATA FILE"))
        self.append(Element('Attribute', name="DATA TYPE"))
        self.append(Element('Attribute', name="NAME"))
