import copy
import itertools
import logging

import cooler
import numpy as np
from matplotlib import cm, colors
from pygenometracks.tracks.GenomeTrack import GenomeTrack
from pygenometracks.utilities import change_chrom_names

DEFAULT_MATRIX_COLORMAP = "RdYlBu_r"
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class HiCMatrixCoolTrack(GenomeTrack):
    SUPPORTED_ENDINGS = [".cool", ".mcool"]
    TRACK_TYPE = "cooler"
    OPTIONS_TXT = f"""
# The different options for color maps can be found here:
# https://matplotlib.org/users/colormaps.html
# the default color map is RdYlBu_r (_r) stands for reverse
# If you want your own colormap you can put the values of the color you want
# For example, colormap = ['blue', 'yellow', 'red']
# or colormap = ['white', (1, 0.88, .66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
#colormap = RdYlBu_r
# depth is the maximum distance that should be plotted.
# If it is more than 125% of the plotted region, it will
# be adjsted to this maximum value.
depth = 100000
# height of track (in cm) can be given.
# Otherwise, the height is computed such that the proportions of the
# hic matrix are kept (e.g. the image does not appear shrink or extended)
# height = 10
# min_value and max_value refer to the contacts in the matrix.
#min_value =2.8
#max_value = 3.0
# the matrix can be transformed using the log1p (or log or -log, but zeros could be problematic)
transform = log1p
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = false
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale_factor = 1
# You can choose to keep the matrix as not rasterized
# (only used if you use pdf or svg output format) by using:
# rasterize = false
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {
        "region": None,  # Cannot be set manually but is set by tracksClass
        "depth": 100000,
        "orientation": None,
        "show_masked_bins": False,
        "scale_factor": 1,
        "transform": "no",
        "max_value": None,
        "min_value": None,
        "rasterize": True,
        "colormap": DEFAULT_MATRIX_COLORMAP,
    }
    NECESSARY_PROPERTIES = ["file"]
    SYNONYMOUS_PROPERTIES = {"max_value": {"auto": None}, "min_value": {"auto": None}}
    POSSIBLE_PROPERTIES = {"orientation": [None, "inverted"], "transform": ["no", "log", "log1p", "-log"]}
    BOOLEAN_PROPERTIES = ["show_masked_bins", "rasterize"]
    STRING_PROPERTIES = ["file", "file_type", "overlay_previous", "orientation", "transform", "title", "colormap"]
    FLOAT_PROPERTIES = {
        "max_value": [-np.inf, np.inf],
        "min_value": [-np.inf, np.inf],
        "scale_factor": [-np.inf, np.inf],
        "height": [0, np.inf],
    }
    INTEGER_PROPERTIES = {"depth": [1, np.inf]}

    # The colormap can only be a colormap

    def __init__(self, *args, **kwargs):
        self.img = None
        self.hic_ma = None
        self.chrom_sizes = None
        self.binsize = None
        self.cmap = None
        self.norm = None

        super().__init__(*args, **kwargs)
        log.debug(f"FILE {self.properties}")

    def set_properties_defaults(self):
        super().set_properties_defaults()
        # Put default img to None for y-axis
        self.img = None
        self.hic_ma = cooler.Cooler(self.properties["file"])
        self.chrom_sizes = self.hic_ma.chromsizes
        self.binsize = self.hic_ma.binsize

        max_depth_in_bins = int(self.properties["depth"] / self.binsize)
        # If the depth is smaller than the binsize. It will display an empty plot
        if max_depth_in_bins < 1:
            self.log.warning(
                f"*Warning*\nThe depth({self.properties['depth']})"
                f" is smaller than binsize({self.binsize})"
                "This will generate an empty track!!\n"
            )
            return

        self.process_color("colormap", colormap_possible=True, colormap_only=True, default_value_is_colormap=True)
        self.cmap = copy.copy(cm.get_cmap(self.properties["colormap"]))
        self.cmap.set_bad("black")

    def plot(self, ax, chrom_region, region_start, region_end):
        log.debug(f"chrom_region {chrom_region}, region_start {region_start}, region_end {region_end}")
        if chrom_region not in self.chrom_sizes:
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.chrom_sizes:
                self.log.warning(
                    f"*Warning*\nNeither {chrom_region_before} nor {chrom_region} exists as a "
                    f"chromosome name on the matrix. This will generate an empty track!!\n"
                )
                self.img = None
                return

        chrom_region = self.check_chrom_str_bytes(self.chrom_sizes, chrom_region)
        if region_end > self.chrom_sizes[chrom_region]:
            self.log.warning(
                "*Warning*\nThe region to plot extends beyond the"
                " chromosome size. Please check.\n"
                f"{chrom_region} size: {self.chrom_sizes[chrom_region]}"
                f". Region to plot {region_start}-{region_end}\n"
            )

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        start_bp = max(0, region_start - self.properties["depth"])
        end_bp = min(self.chrom_sizes[chrom_region], region_end + self.properties["depth"])
        p1 = int(np.floor(start_bp / self.binsize))
        p2 = int(np.ceil(end_bp / self.binsize))

        # start pos need to include the last end
        start_pos = [int(i * self.binsize) for i in range(p1, p2 + 1)]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix(balance=False, sparse=False).fetch(f"{chrom_region}:{start_bp}-{end_bp}")
        matrix = np.triu(matrix)

        # limit the 'depth' based on the length of the region being viewed
        region_len = region_end - region_start
        depth = min(self.properties["depth"], int(region_len * 1.25))
        if depth < self.properties["depth"]:
            log.warning(
                f"The depth was set to {self.properties['depth']} which is more than 125%"
                " of the region plotted. The depth will be set "
                f"to {depth}.\n"
            )

        # Replace all nan values by 0.
        np.nan_to_num(matrix, nan=0, copy=False)

        # scale
        matrix = matrix * self.properties["scale_factor"]

        if self.properties["transform"] == "log1p":
            matrix += 1

        elif self.properties["transform"] in ["-log", "log"]:
            if matrix.min() < 0:
                raise ValueError("HiC Matrix contains negative value, can not perform log transform.")

            # We first replace 0 values by minimum values after 0
            mask = matrix == 0
            try:
                matrix[mask] = matrix[~mask].min()
                matrix = np.log(matrix)
            except ValueError:
                self.log.info("All values are 0, no log applied.")
            else:
                if self.properties["transform"] == "-log":
                    matrix *= -1

        if self.properties["max_value"] is not None:
            vmax = self.properties["max_value"]
        else:
            # try to use a 'aesthetically pleasant' max value
            try:
                vmax = np.percentile(matrix.diagonal(1), 80)
            except IndexError:
                vmax = None

        if self.properties["min_value"] is not None:
            vmin = self.properties["min_value"]
        else:
            vmin = matrix.min()
            # if the region length is large with respect to the chromosome length, the diagonal may have
            # very few values or none. Thus, the following lines reduce the number of bins until the
            # diagonal is at least length 5 but make sure you have at least one value:
            num_bins_from_diagonal = max(1, int(region_len / self.binsize))
            for num_bins in range(0, num_bins_from_diagonal)[::-1]:
                distant_diagonal_values = matrix.diagonal(num_bins)
                if len(distant_diagonal_values) > 5:
                    break
                vmin = min(vmin, np.median(distant_diagonal_values))
        self.log.info(
            "setting min, max values for track " f"{self.properties['section_name']} to: " f"{vmin}, {vmax}\n"
        )
        if self.properties["transform"] == "log1p":
            self.norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            self.norm = colors.Normalize(vmin=vmin, vmax=vmax)

        self.img = self.pcolormesh_45deg(ax, matrix, start_pos)
        if self.properties["rasterize"]:
            self.img.set_rasterized(True)
        if self.properties["orientation"] == "inverted":
            ax.set_ylim(depth, 0)
        else:
            ax.set_ylim(0, depth)

    def plot_y_axis(self, cbar_ax, plot_ax):
        if self.img is None:
            return

        GenomeTrack.plot_custom_cobar(self, cbar_ax)

    def pcolormesh_45deg(self, ax, matrix_c, start_pos_vector):
        """
        Plot a 45 degree heatmap.

        Turns the matrix 45 degrees and adjusts the bins to match the actual start end positions.
        """
        # code for rotating the image 45 degrees
        n = matrix_c.shape[0]
        # create rotation/scaling matrix
        t = np.array([[1, 0.5], [-1, 0.5]])
        # create coordinate matrix and transform it
        matrix_a = np.dot(
            np.array([(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)]), t
        )
        # this is to convert the indices into bp ranges
        x = matrix_a[:, 1].reshape(n + 1, n + 1)
        y = matrix_a[:, 0].reshape(n + 1, n + 1)
        # plot
        im = ax.pcolormesh(x, y, np.flipud(matrix_c), cmap=self.cmap, norm=self.norm)
        return im
