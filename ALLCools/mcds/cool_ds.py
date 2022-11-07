import json
import pathlib

import numpy as np
import pandas as pd
import xarray as xr
from cooler import read_chromsizes
from scipy import ndimage


class CoolDS:
    def __init__(self, cool_ds_paths, chrom_sizes_path, sample_weights: pd.Series, sample_dim="sample_id"):
        """
        Multiple chromatin conformation matrix profiles

        Parameters
        ----------
        cool_ds_paths
            List of paths to cool ds or path wildcard
        chrom_sizes_path
            Path to chrom sizes file
        sample_weights: pd.Series
            A series of sample weights used for sum sample
        sample_dim
            Name of sample dimension
        """
        self.cool_ds_paths = self._prepare_cool_ds_paths(cool_ds_paths)
        self.chrom_sizes_path = chrom_sizes_path
        self.chrom_sizes = read_chromsizes(self.chrom_sizes_path)
        self.bin_size = self._get_bin_size_from_cool_ds()
        self.sample_weights = sample_weights
        self.sample_dim = sample_dim
        self.sample_weights.index.name = self.sample_dim

        self._chrom_ds_cache = {}

        return

    def _get_bin_size_from_cool_ds(self):
        _path = self.cool_ds_paths[0]
        _chrom = self.chrom_sizes.index[0]
        with open(f"{_path}/{_chrom}/.zattrs") as f:
            attrs = json.load(f)
        return attrs["cooler_bin_size"]

    @staticmethod
    def _prepare_cool_ds_paths(cool_ds_paths):
        if isinstance(cool_ds_paths, (str, pathlib.Path)):
            if "*" in str(cool_ds_paths):

                import glob

                cool_ds_paths = sorted(glob.glob(cool_ds_paths))
            else:
                cool_ds_paths = [cool_ds_paths]
        return cool_ds_paths

    def load_chrom_ds(self, chrom, chrom2=None):
        """
        Load chrom ds from cache or disk

        Parameters
        ----------
        chrom
            Chromosome name
        chrom2
            Chromosome name for second dimension, optional

        Returns
        -------
        xarray.Dataset
        """
        if (chrom, chrom2) in self._chrom_ds_cache:
            return self._chrom_ds_cache[(chrom, chrom2)]
        else:
            if chrom not in self.chrom_sizes.index:
                raise ValueError(f"Chrom {chrom} not found in chrom_sizes")

            if chrom2 is None:
                chrom_paths = [f"{path}/{chrom}" for path in self.cool_ds_paths]
            else:
                raise NotImplementedError

            chrom_ds = xr.open_mfdataset(
                chrom_paths,
                concat_dim="sample_id",
                combine="nested",
                engine="zarr",
            )
            self._chrom_ds_cache[(chrom, chrom2)] = chrom_ds
            return chrom_ds

    def _region_to_slice(self, chrom, start=None, end=None):
        assert chrom in self.chrom_sizes.index, f"Chrom {chrom} not found in chrom_sizes"
        if start is not None:
            start_bin = start // self.bin_size
        else:
            start_bin = 0
        if end is not None:
            end_bin = end // self.bin_size + 1
        else:
            end_bin = self.chrom_sizes[chrom] // self.bin_size + 1
        return slice(start_bin, end_bin)

    def fetch(self, chrom, start=None, end=None, chrom2=None, start2=None, end2=None):
        """
        Fetch single chromosome matrix dataset from cool ds

        Parameters
        ----------
        chrom
            Chromosome name
        start
            Start position
        end
            End position
        chrom2
            Chromosome name for second dimension, optional
        start2
            Start position for second dimension, optional
        end2
            End position for second dimension, optional

        Returns
        -------
        CoolDSChrom
        """
        _chrom_ds = self.load_chrom_ds(chrom, chrom2)
        if chrom2 is None:
            chrom2 = chrom

        slice1 = self._region_to_slice(chrom, start, end)
        if start2 is None and end2 is None:
            slice2 = slice1
        else:
            slice2 = self._region_to_slice(chrom2, start2, end2)
        use_ds = _chrom_ds.isel(bin1=slice1, bin2=slice2)
        use_ds.attrs["chrom1_offset"] = 0 if slice1.start is None else slice1.start
        use_ds.attrs["chrom2_offset"] = 0 if slice2.start is None else slice2.start
        return CoolDSChrom(
            use_ds, sample_dim=self.sample_dim, sample_weights=self.sample_weights, bin_size=self.bin_size
        )


class CoolDSChrom(xr.Dataset):
    __slots__ = ()

    def __init__(self, dataset, sample_dim, sample_weights, bin_size):
        super().__init__(data_vars=dataset.data_vars, coords=dataset.coords, attrs=dataset.attrs)
        self.attrs["sample_dim"] = sample_dim

        if bin_size is not None:
            self.attrs["bin_size"] = bin_size

        sample_weights.index.name = sample_dim
        self.coords["sample_weights"] = sample_weights
        return

    @property
    def bin_size(self):
        return self.attrs["bin_size"]

    @property
    def sample_dim(self):
        return self.attrs["sample_dim"]

    @property
    def sample_weights(self):
        if "sample_weights" in self.coords:
            return self.coords["sample_weights"]
        else:
            raise KeyError("sample_weights not found in coords")

    def matrix(
        self,
        samples,
        value_type,
        scale_factor=10e6,
        fill_lower_triangle=False,
        log1p=True,
        da_name="real",
        rotate=False,
        rotate_cval=np.NaN,
        rotate_height_bp=5000000,
        dtype="float32",
    ):
        """
        Fetch matrix from cool ds

        Parameters
        ----------
        samples
            Samples to fetch, if multiple samples are provided, the matrix is summed with sample weights
        value_type
            Value type to fetch
        scale_factor
            Scale factor to apply to matrix
        fill_lower_triangle
            Fill matrix lower triangle with upper triangle
        log1p
            Apply log1p to matrix
        da_name
            DataArray name
        rotate
            Rotate matrix by 45 degrees to make triangle plot
        rotate_cval
            Rotate matrix fill value
        rotate_height_bp
            Rotate matrix height in base pairs
        dtype
            Data type of the returned matrix

        Returns
        -------
        np.ndarray
        """
        sample_da = self.sel({self.sample_dim: samples, f"{da_name}_value_type": value_type})[da_name]
        if not isinstance(samples, str):
            # sum sample_dim if multiple samples selected
            sample_da = (sample_da * self.sample_weights).sum(dim=self.sample_dim) / (
                self.sample_weights.sum() / scale_factor
            )

        data = sample_da.values.astype(dtype)

        if fill_lower_triangle:
            # complete the lower triangle
            data = data + data.T - np.diag(data.diagonal())

        if log1p:
            data = np.log1p(data)

        if rotate:
            # to make triangle plot, rotate the matrix by 45 degrees
            data = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False, cval=rotate_cval)

            middle = data.shape[0] // 2
            height = (rotate_height_bp // self.bin_size) / np.sqrt(2) + 4
            height = int(height)
            height = min(height, data.shape[0] // 2)

            bottom = middle - height
            if fill_lower_triangle:
                top = middle + height
            else:
                top = middle
            data = data[bottom:top].copy()
        return data
