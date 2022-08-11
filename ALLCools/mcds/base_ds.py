import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from ALLCools.utilities import parse_mc_pattern


def _chunk_pos_to_bed_df(chrom, chunk_pos):
    records = []
    for i in range(len(chunk_pos) - 1):
        start = chunk_pos[i]
        end = chunk_pos[i + 1]
        records.append([chrom, start, end])
    return pd.DataFrame(records, columns=["chrom", "start", "end"])


class Codebook(xr.DataArray):
    """The Codebook data array records methyl-cytosine context in genome."""

    __slots__ = ()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # in-memory cache for the positions of cytosines matching the pattern
        self.attrs["__mc_pos_cache"] = {}
        self.attrs["__mc_pos_bool_cache"] = {}

        # the continuity of codebook is determined by occurrence of pos coords
        # and should not be changed
        self.attrs["__continuity"] = False if "pos" in self.coords else True

    @property
    def continuity(self):
        """Whether the codebook is continuous or not."""
        return self.attrs["__continuity"]

    @property
    def mc_type(self):
        """The type of methyl-cytosine context."""
        return self.get_index("mc_type")

    def get_mc_pos_bool(self, mc_pattern):
        """Get the boolean array of cytosines matching the pattern."""
        if mc_pattern in self.attrs["__mc_pos_bool_cache"]:
            return self.attrs["__mc_pos_bool_cache"][mc_pattern]
        else:
            if mc_pattern is None:
                # get all cytosines
                judge = np.ones_like(self.mc_type, dtype=bool)
            else:
                # get cytosines matching the pattern
                judge = self.mc_type.isin(parse_mc_pattern(mc_pattern))
            _bool = self.sel(mc_type=judge).sum(dim="mc_type").values.astype(bool)
            self.attrs["__mc_pos_bool_cache"][mc_pattern] = _bool
            return _bool

    def get_mc_pos(self, mc_pattern):
        """Get the positions of cytosines matching the pattern."""
        if mc_pattern in self.attrs["__mc_pos_cache"]:
            return self.attrs["__mc_pos_cache"][mc_pattern]
        else:
            _bool = self.get_mc_pos_bool(mc_pattern)

            if self.continuity:
                _pos = self.coords["pos"].values[_bool]
            else:
                _pos = np.where(_bool)[0]

            self.attrs["__mc_pos_cache"][mc_pattern] = _pos
            return _pos


class BaseDSChrom(xr.Dataset):
    """The BaseDS class for data within single chromosome."""

    __slots__ = ()

    def __init__(self, dataset, coords=None, attrs=None):
        if isinstance(dataset, xr.Dataset):
            data_vars = dataset.data_vars
            coords = dataset.coords if coords is None else coords
            attrs = dataset.attrs if attrs is None else attrs
        else:
            data_vars = dataset
        super().__init__(data_vars=data_vars, coords=coords, attrs=attrs)

        # continuity is set to True by default
        self.continuity = True
        return

    @property
    def continuity(self):
        """
        The continuity of pos dimension.

        the BaseDSChrom has tow mode on the position dimension:
        1. "continuous" mode: the position dimension is continuous
        and all bases (including non-cytosines) are recorded.
        In this mode, the position dimension do not have coordinates, the index is genome position.
        2. "discrete" mode: the position dimension is discrete due to some discontinuous selection.
        In this mode, the position dimension has coordinates.
        """
        return self.attrs["__continuity"]

    @continuity.setter
    def continuity(self, value):
        """Set the continuity of pos dimension."""
        if not value:
            assert "pos" in self.coords, (
                "The position dimension is set to discontinuous, " "but the pos coords is missing."
            )
            self.attrs["__continuity"] = False
        else:
            self.attrs["__continuity"] = True

    def _clear_attr_cache(self):
        """Clear the attr cache."""
        for attr in list(self.attrs.keys()):
            if str(attr).startswith("__"):
                del self.attrs[attr]

    def _continuous_pos_selection(self, start, end):
        """Select the positions to create a continuous BaseDSChrom."""
        # once the pos is selected, the ds is continuous
        # one must set the pos coords and set the continuity to True
        if start is not None or end is not None:
            obj = self.sel(pos=slice(start, end, None))
            return obj
        else:
            return self

    def _discontinuous_pos_selection(self, pos_sel):
        """Select the positions to create a discontinuous BaseDSChrom."""
        # once the pos is selected, the ds is not continuous anymore
        # one must set the pos coords and set the continuity to False
        self._clear_attr_cache()
        ds = self.sel(pos=pos_sel).assign_coords(pos=pos_sel)
        ds.continuity = False
        return ds

    @classmethod
    def open(cls, path, start=None, end=None):
        """
        Open a BaseDSChrom object from a zarr path.

        If start and end are not None, only the specified region will be opened.

        Parameters
        ----------
        path
            The zarr path to the chrom dataset.
        start
            The start position of the region to be opened.
        end
            The end position of the region to be opened.

        Returns
        -------
        BaseDSChrom
        """
        path = pathlib.Path(path)
        _obj = cls(xr.open_zarr(path, decode_cf=False))
        _obj = _obj._continuous_pos_selection(start, end)
        return _obj

    @property
    def chrom(self):
        """The chromosome name."""
        return self.attrs["chrom"]

    @property
    def chrom_size(self):
        """The chromosome size."""
        return self.attrs["chrom_size"]

    @property
    def obs_dim(self):
        """The observation dimension name."""
        return self.attrs["obs_dim"]

    @property
    def obs_size(self):
        """The observation size."""
        return self.attrs["obs_size"]

    @property
    def obs_names(self):
        """The observation names."""
        return self.get_index(self.obs_dim)

    @property
    def mc_types(self):
        """The methyl-cytosine types."""
        return self.get_index("mc_type")

    @property
    def chrom_chunk_pos(self):
        """The chromosome chunk position."""
        return self.get_index("chunk_pos")

    @property
    def chrom_chunk_bed_df(self) -> pd.DataFrame:
        """The chromosome chunk bed dataframe."""
        chunk_pos = self.chrom_chunk_pos
        bed = _chunk_pos_to_bed_df(self.chrom, chunk_pos)
        return bed

    @property
    def codebook(self) -> Codebook:
        """Get the codebook data array."""
        # catch the codebook in the attrs, only effective in memory
        if "__cb_obj" not in self.attrs:
            self.attrs["__cb_obj"] = Codebook(self["codebook"])
        return self.attrs["__cb_obj"]

    @property
    def cb(self) -> Codebook:
        """Alias for codebook."""
        return self.codebook

    def select_mc_type(self, pattern):
        cb = self.codebook
        pattern_bool = cb.get_mc_pos(pattern)

        ds = self._discontinuous_pos_selection(pattern_bool)
        return ds
