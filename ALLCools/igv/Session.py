from .BaseClass import *
from .Track import *


def get_panel_int(panel):
    if panel.lower().startswith('d'):
        return 1
    elif panel.lower().startswith('f'):
        return 2
    else:
        raise ValueError('Panel can only be "Data" or "Feature"')


class Session(Element):
    def __init__(self, genome="mm10",
                 hasGeneTrack="false",
                 hasSequenceTrack="true",
                 locus="chr1:19000000-19100000",
                 version="8", refseq_track=True):
        _tag = 'Session'
        _attrib = {'genome': genome,
                   'hasGeneTrack': hasGeneTrack,
                   'hasSequenceTrack': hasSequenceTrack,
                   'locus': locus,
                   'version': version}
        super().__init__(_tag, attrib=_attrib)

        # add children
        self.__add_resources()
        self.__add_data_panel()
        self.__add_feature_panel()
        self.__add_other_routine()

        # add grandchildren
        self.__add_sequence_track()
        if refseq_track:
            self.__add_refseq_gene()
        return

    def __add_resources(self):
        self.append(Resources())

    def __add_data_panel(self):
        self.append(Panel(height="2080", name="DataPanel", width="1903"))

    def __add_feature_panel(self):
        self.append(Panel(height="95", name="FeaturePanel", width="1903"))

    def __add_other_routine(self):
        self.append(PanelLayout())
        self.append(HiddenAttributes())

    def __add_resource(self, path):
        # self[0] is Resources
        self[0].append(Resource(path))

    def __add_sequence_track(self):
        # self[2] is FeaturePanel
        self[2].append(Track(clazz="org.broad.igv.track.SequenceTrack",
                             fontSize="10",
                             id="Reference sequence",
                             name="Reference sequence",
                             visible="true"))

    def __add_refseq_gene(self):
        self[2].append(Track(altColor="0,0,178",
                             clazz="org.broad.igv.track.FeatureTrack",
                             color="0,0,178",
                             colorScale="ContinuousColorScale;0.0;185.0;255,255,255;0,0,178",
                             fontSize="10",
                             height="35",
                             id="mm10_genes",
                             name="Refseq genes",
                             visible="true"))

    def add_data_track(self, track_kws, data_range_kws, panel='Data'):
        _panel = get_panel_int(panel)
        track = DataTrack(**track_kws)
        track.add_data_range(**data_range_kws)

        self[_panel].append(track)
        # id is track path, checked key existence in Track.__init__
        self.__add_resource(track_kws['id'])
        pass

    def add_feature_track(self, track_kws, panel='Feature'):
        _panel = get_panel_int(panel)
        track = FeatureTrack(**track_kws)
        self[_panel].append(track)
        self.__add_resource(track_kws['id'])
        pass
