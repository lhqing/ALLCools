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

    @property
    def c_pos(self):
        """The positions of cytosines in mC type."""
        return self.attrs["c_pos"]

    @property
    def context_size(self):
        """The size of context in mC type."""
        return self.attrs["context_size"]

    def _validate_mc_pattern(self, mc_pattern):
        mc_pattern = mc_pattern.upper()

        if len(mc_pattern) != self.context_size:
            raise ValueError(
                f"The length of mc_pattern {len(mc_pattern)} is not equal to context_size {self.context_size}."
            )
        if mc_pattern[self.c_pos] != "C":
            raise ValueError(f"The c_pos position {self.c_pos} in mc_pattern {mc_pattern} is not a cytosine.")

        return mc_pattern

    def get_mc_pos_bool(self, mc_pattern):
        """Get the boolean array of cytosines matching the pattern."""
        mc_pattern = self._validate_mc_pattern(mc_pattern)

        if mc_pattern in self.attrs["__mc_pos_bool_cache"]:
            return self.attrs["__mc_pos_bool_cache"][mc_pattern]
        else:
            if mc_pattern is None:
                # get all mc types
                judge = np.ones_like(self.mc_type, dtype=bool)
            else:
                # get mc types matching the pattern
                judge = self.mc_type.isin(parse_mc_pattern(mc_pattern))
            _bool = self.sel(mc_type=judge).sum(dim="mc_type").values.astype(bool)
            self.attrs["__mc_pos_bool_cache"][mc_pattern] = _bool
            return _bool

    def get_mc_pos(self, mc_pattern):
        """Get the positions of mc types matching the pattern."""
        mc_pattern = self._validate_mc_pattern(mc_pattern)

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
        self.offset = 0
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
            # when continuity is set to False, the offset is set to None
            self.offset = None
            self.attrs["__continuity"] = False
        else:
            self.attrs["__continuity"] = True

    @property
    def offset(self):
        """The offset of the position dimension, only valid when continuity is True."""
        return self.attrs["__offset"]

    @offset.setter
    def offset(self, value):
        """Set the offset of the position dimension."""
        if value is not None and not self.continuity:
            raise ValueError("The offset is only valid when the position dimension is continuous.")
        self.attrs["__offset"] = value

    def clear_attr_cache(self):
        """Clear the attr cache."""
        for attr in list(self.attrs.keys()):
            if str(attr).startswith("__"):
                del self.attrs[attr]

    def _continuous_pos_selection(self, start, end):
        """Select the positions to create a continuous BaseDSChrom."""
        # for continuous mode, the pos should have an offset to convert to genome position
        if not self.continuity:
            raise ValueError("The position dimension is not continuous, unable to perform _continuous_pos_selection.")

        if start is not None or end is not None:
            if start is not None:
                start -= self.offset
            if end is not None:
                end -= self.offset
            obj = self.sel(pos=slice(start, end, None))
            if start is not None:
                obj.offset = start
            return obj
        else:
            return self

    def _discontinuous_pos_selection(self, pos_sel=None, idx_sel=None):
        """
        Select the positions to create a discontinuous BaseDSChrom.

        Parameters
        ----------
        pos_sel :
            using genome position to select the positions,
            for continuous mode, the pos should have an offset to convert to idx position
        idx_sel :
            using idx position to select the positions

        Returns
        -------
        BaseDSChrom
        """
        # once the pos is selected, the ds is not continuous anymore
        # one must set the pos coords and set the continuity to False
        if idx_sel is not None and pos_sel is not None:
            raise ValueError("Only one of idx_sel and pos_sel can be specified.")
        elif idx_sel is not None:
            pass
        elif pos_sel is not None:
            if self.continuity:
                # da is continuous, convert pos to idx
                idx_sel = pos_sel - self.offset
            else:
                # da is not continuous, treat pos_sel as idx_sel
                idx_sel = pos_sel
        else:
            raise ValueError("One of idx_sel or pos_sel must be specified.")

        if self.continuity:
            # if the ds was continuous, add offset to the pos coords and turn off continuity
            offset_to_add = self.offset
        else:
            offset_to_add = 0

        ds = self.sel(pos=idx_sel).assign_coords(pos=idx_sel + offset_to_add)
        ds.clear_attr_cache()
        ds.continuity = False
        return ds

    @classmethod
    def open(cls, path, start=None, end=None, codebook_path=None):
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
        codebook_path
            The path to the codebook file if the BaseDS does not have a codebook.
            Codebook contexts, c_pos, and shape must be compatible with the BaseDS.

        Returns
        -------
        BaseDSChrom
        """
        path = pathlib.Path(path)

        _zarr_obj = xr.open_zarr(path, decode_cf=False)

        if "codebook" not in _zarr_obj.data_vars:
            if codebook_path is None:
                raise ValueError("The BaseDS does not have a codebook, but no codebook_path is specified.")
            _cb = xr.open_zarr(codebook_path, decode_cf=False)["codebook"]
            # validate _cb attrs compatibility
            flag = True
            _cb_mc_types = _cb.get_index("mc_type").values
            _obj_mc_types = _zarr_obj.get_index("mc_type").values
            # noinspection PyUnresolvedReferences
            _diff = (_cb_mc_types != _obj_mc_types).sum()
            if _diff > 0:
                flag = False
                print("The codebook mc_types are not compatible with the BaseDS.")
            if _cb.shape[0] != _zarr_obj["data"].shape[0]:
                flag = False
                print("The codebook shape is not compatible with the BaseDS.")
            if not flag:
                raise ValueError("The BaseDS and codebook are not compatible.")
            _zarr_obj["codebook"] = _cb

        _obj = cls(_zarr_obj)
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
            self.attrs["__cb_obj"].attrs["c_pos"] = self.attrs["c_pos"]
            self.attrs["__cb_obj"].attrs["context_size"] = self.attrs["context_size"]
        return self.attrs["__cb_obj"]

    @property
    def cb(self) -> Codebook:
        """Alias for codebook."""
        return self.codebook

    def select_mc_type(self, pattern):
        cb = self.codebook
        pattern_bool = cb.get_mc_pos(pattern)

        ds = self._discontinuous_pos_selection(idx_sel=pattern_bool)
        return ds

    @property
    def pos_index(self):
        """The position index."""
        return self.get_index("pos")
