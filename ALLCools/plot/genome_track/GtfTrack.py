"""
This GtfTrack class replaced the original pygenometracks class.

It supports read-in a gffutils database, rather than parse very time
"""

import collections
import warnings

import gffutils
import numpy as np
from matplotlib import font_manager
from pygenometracks.tracks.BedTrack import BedTrack
from pygenometracks.tracks.GenomeTrack import GenomeTrack
from pygenometracks.utilities import InputError, temp_file_from_intersect

DEFAULT_BED_COLOR = "#1f78b4"
DISPLAY_BED_VALID = ["collapsed", "triangles", "interleaved", "stacked"]
DISPLAY_BED_SYNONYMOUS = {"interlaced": "interleaved", "domain": "interleaved"}
DEFAULT_DISPLAY_BED = "stacked"
AROUND_REGION = 100000


def _is_sqlite3(path):
    with open(path, "rb") as f:
        head = f.read(16)
        if head.decode("utf8").startswith("SQLite format"):
            return True
        else:
            return False


warnings.filterwarnings(
    "ignore",
    message="It appears you have a gene feature"
    " in your GTF file. You may want to use the "
    "`disable_infer_genes` option to speed up database "
    "creation",
)
warnings.filterwarnings(
    "ignore",
    message="It appears you have a transcript "
    "feature in your GTF file. You may want to use the "
    "`disable_infer_transcripts` option to speed up "
    "database creation",
)
# In gffutils v0.10 they changed the error message:
warnings.filterwarnings(
    "ignore",
    message="It appears you have a gene feature"
    " in your GTF file. You may want to use the "
    "`disable_infer_genes=True` option to speed up database "
    "creation",
)
warnings.filterwarnings(
    "ignore",
    message="It appears you have a transcript "
    "feature in your GTF file. You may want to use the "
    "`disable_infer_transcripts=True` option to speed up "
    "database creation",
)


class ReadGtf:
    """Read a gtf file."""

    def __init__(self, file_path, prefered_name="transcript_name", merge_transcripts=True):
        self.file_type = "bed12"

        # list of bed fields
        self.fields = [
            "chromosome",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "thick_start",
            "thick_end",
            "rgb",
            "block_count",
            "block_sizes",
            "block_starts",
        ]

        self.BedInterval = collections.namedtuple("BedInterval", self.fields)
        # I think the name which should be written
        # should be the transcript_name
        # But we can change it to gene_name
        self.prefered_name = prefered_name
        self.merge_transcripts = merge_transcripts

        # Will process the gtf to get one item per transcript:
        # This will create a database:
        try:
            if _is_sqlite3(file_path):
                self.db = gffutils.FeatureDB(file_path)
            else:
                self.db = gffutils.create_db(file_path, ":memory:")
        except ValueError as ve:
            if "No lines parsed" in str(ve):
                self.length = 0
                self.all_transcripts = open(file_path)
            else:
                raise InputError("This is not a gtf file.")
        else:
            if self.merge_transcripts:
                self.length = len(list(self.db.features_of_type("gene")))
                self.all_transcripts = self.db.features_of_type("gene", order_by="start")
            else:
                self.length = len(list(self.db.features_of_type("transcript")))
                self.all_transcripts = self.db.features_of_type("transcript", order_by="start")

    def __iter__(self):
        return self

    def __next__(self):
        """Get next bedInterval object."""
        bed = self.get_bed_interval()

        return bed

    def get_bed_interval(self):
        """
        Get the bed interval.

        Process a transcript from the database, retrieve all the values and return a namedtuple object.
        """
        tr = next(self.all_transcripts)
        # The name would be the prefered_name if exists
        try:
            trName = tr.attributes[self.prefered_name][0]
        except KeyError:
            # Else try to guess the prefered_name from exons:
            try:
                trName = {
                    e.attributes[self.prefered_name][0]
                    for e in self.db.children(tr, featuretype="exon", order_by="start")
                }.pop()
            except KeyError:
                # Else take the transcript id
                trName = tr.id
        # If the cds is defined in the gtf,
        # use it to define the thick start and end
        # The gtf is 1-based closed intervalls
        # and bed are 0-based half-open so:
        # I need to remove one from each start
        try:
            cds_start = next(self.db.children(tr, featuretype="CDS", order_by="start")).start - 1
            cds_end = next(self.db.children(tr, featuretype="CDS", order_by="-start")).end
        except StopIteration:
            # If the CDS is not defined, then it is set to the start
            # as proposed here:
            # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            cds_start = tr.start - 1
            cds_end = tr.start - 1
        # Get all exons starts and end to get lengths
        exons_starts = [e.start - 1 for e in self.db.children(tr, featuretype="exon", order_by="start")]
        exons_ends = [e.end for e in self.db.children(tr, featuretype="exon", order_by="start")]
        exons_length = [e - s for s, e in zip(exons_starts, exons_ends)]
        relative_exons_starts = [s - (tr.start - 1) for s in exons_starts]
        line_values = [
            tr.chrom,
            tr.start - 1,
            tr.end,
            trName,
            0,
            tr.strand,
            cds_start,
            cds_end,
            "0",
            len(exons_starts),
            exons_length,
            relative_exons_starts,
        ]
        return self.BedInterval._make(line_values)


class GtfTrack(BedTrack):
    """GTF track."""

    SUPPORTED_ENDINGS = ["gtf", "gtf.gz", "gtf.db"]
    TRACK_TYPE = "gtf"
    OPTIONS_TXT = (
        GenomeTrack.OPTIONS_TXT
        + f"""
# By default the transcript_name is used.
# If you want to use the gene_name:
# prefered_name = gene_name
# By default, the gtf is transformed to transcripts
# If you want to use see only one structure per gene
# merge_transcripts = true
# You can change the color of coding sequences by:
color = darkblue
# height of track in cm
height = 5
# whether printing the labels
labels = false
# optional:
# by default the labels are not printed if you have more than 60 features.
# to change it, just increase the value:
#max_labels = 60
# optional: font size can be given to override the default size
fontsize = 10
# optional: line_width
#line_width = 0.5
# the display parameter defines how the gtf file is plotted.
# Default is 'stacked' where regions are plotted on different lines so
# we can see all regions and all labels.
# The other options are ['collapsed', 'interleaved', 'triangles']
# These options assume that the regions do not overlap.
# `collapsed`: The gtf regions are plotted one after the other in one line.
# `interleaved`: The gtf regions are plotted in two lines, first up, then down, then up etc.
# optional, default is black. To remove the border, simply set 'border_color' to none
# Not used in tssarrow style
#border_color = black
# style to plot the genes when the display is not triangles
#style = UCSC
#style = flybase
#style = tssarrow
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size
# otherwise, on each image the height of each gene changes
#gene_rows = 10
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global_max_row = true
# If you want to plot all labels inside the plotting region:
#all_labels_inside = true
# If you want to display the name of the gene which goes over the plotted
# region in the right margin put:
#labels_in_margin = true
# if you use UCSC style, you can set the relative distance between 2 arrows on introns
# default is 2
#arrow_interval = 2
# if you use tssarrow style, you can choose the length of the arrow in bp
# (default is 4% of the plotted region)
#arrow_length = 5000
# if you use flybase or tssarrow style, you can choose the color of non-coding intervals:
#color_utr = grey
# as well as the proportion between their height and the one of coding
# (by default they are the same height):
#height_utr = 1
# By default, for oriented intervals in flybase style,
# or bed files with less than 12 columns, the arrowhead is added
# outside of the interval.
# If you want that the tip of the arrow correspond to
# the extremity of the interval use:
# arrowhead_included = true
# optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """
    )

    DEFAULTS_PROPERTIES = {
        "fontsize": 12,
        "orientation": None,
        "color": DEFAULT_BED_COLOR,
        "border_color": "black",
        "labels": True,
        "style": "flybase",
        "display": DEFAULT_DISPLAY_BED,
        "line_width": 0.5,
        "max_labels": 60,
        "prefered_name": "transcript_name",
        "merge_transcripts": False,
        "global_max_row": False,
        "gene_rows": None,
        "arrow_interval": 2,
        "arrowhead_included": False,
        "color_utr": "grey",
        "height_utr": 1,
        "arrow_length": None,
        "region": None,  # Cannot be set manually but is set by tracksClass
        "all_labels_inside": False,
        "labels_in_margin": False,
    }
    NECESSARY_PROPERTIES = ["file"]
    SYNONYMOUS_PROPERTIES = {"display": DISPLAY_BED_SYNONYMOUS}
    POSSIBLE_PROPERTIES = {
        "orientation": [None, "inverted"],
        "style": ["flybase", "UCSC", "tssarrow"],
        "display": DISPLAY_BED_VALID,
    }
    BOOLEAN_PROPERTIES = [
        "labels",
        "merge_transcripts",
        "global_max_row",
        "arrowhead_included",
        "all_labels_inside",
        "labels_in_margin",
    ]
    STRING_PROPERTIES = [
        "prefered_name",
        "file",
        "file_type",
        "overlay_previous",
        "orientation",
        "title",
        "style",
        "color",
        "border_color",
        "color_utr",
        "display",
    ]
    FLOAT_PROPERTIES = {"fontsize": [0, np.inf], "line_width": [0, np.inf], "height": [0, np.inf], "height_utr": [0, 1]}
    INTEGER_PROPERTIES = {
        "gene_rows": [0, np.inf],
        "max_labels": [0, np.inf],
        "arrow_interval": [1, np.inf],
        "arrow_length": [0, np.inf],
    }

    def set_properties_defaults(self):
        """Set the default values for the properties."""
        super(BedTrack, self).set_properties_defaults()
        self.fp = font_manager.FontProperties(size=self.properties["fontsize"])
        self.colormap = None
        # check if the color given is a color map
        # Contrary to bed it cannot be a colormap
        self.process_color("color", colormap_possible=False, bed_rgb_possible=False, default_value_is_colormap=False)

        # check if border_color and color_utr are colors
        # if they are part of self.properties
        # (for example, TADsTracks do not have color_utr)
        for param in [p for p in ["border_color", "color_utr"] if p in self.properties]:
            self.process_color(param, bed_rgb_possible=False)

        # to set the distance between rows
        self.row_scale = 2.3

    def get_bed_handler(self, plot_regions=None):
        """Get the bed handler for the track."""
        file_to_open = self.properties["file"]
        if not self.properties["global_max_row"]:
            if plot_regions is not None:
                if _is_sqlite3(self.properties["file"]):
                    print("global_max_row not supported when gtf is provided as SQLite3 db")
                else:
                    # I do the intersection:
                    file_to_open = temp_file_from_intersect(self.properties["file"], plot_regions, AROUND_REGION)

        gtf_db = ReadGtf(file_to_open, self.properties["prefered_name"], self.properties["merge_transcripts"])
        total_length = gtf_db.length
        return gtf_db, total_length
