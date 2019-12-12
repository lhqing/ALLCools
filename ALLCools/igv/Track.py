from .BaseClass import Track, DataRange
from .defaults import *


class DataTrack(Track):
    def __init__(self, **kwargs):
        clazz = 'org.broad.igv.track.DataSourceTrack'
        kwargs['clazz'] = clazz

        _attrib = DATA_TRACK_DEFAULT.copy()
        _attrib.update(kwargs)
        super().__init__(**_attrib)

        self.add_data_range()  # add a default data range

    def add_data_range(self, **kwargs):
        # this is only one data range allowed
        for data_range in self.findall('DataRange'):
            self.remove(data_range)
        self.append(DataRange(**kwargs))


class FeatureTrack(Track):
    def __init__(self, **kwargs):
        clazz = 'org.broad.igv.track.FeatureTrack'
        kwargs['clazz'] = clazz

        _attrib = FEATURE_TRACK_DEFAULT.copy()
        _attrib.update(kwargs)
        super().__init__(**_attrib)
