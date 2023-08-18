import json
import pathlib
from typing import Union

import dask
import numpy as np
import pandas as pd
import xarray as xr
from cooler import binnify, read_chromsizes
from scipy import ndimage


def _get_chrom_offsets(bins):
    _co = {chrom: bins[bins["chrom"] == chrom].index[0] for chrom in bins["chrom"].cat.categories}
    return _co


class CoolDS:
    def __init__(
        self,
        cool_ds_paths,
        chrom_sizes_path,
        sample_weights: Union[pd.Series, None] = None,
        sample_dim="sample_id",
        bin_size=None,
    ):
        """
        Multiple chromatin conformation matrix profiles.

        Parameters
        ----------
        cool_ds_paths
            List of paths to cool ds or path wildcard
        chrom_sizes_path
            Path to chrom sizes file
        sample_weights: pd.Series, None
            A series of sample weights used for sum sample
        sample_dim
            Name of sample dimension
        bin_size
            Bin size, if None, will be inferred from cool ds
        """
        self.cool_ds_paths = self._prepare_cool_ds_paths(cool_ds_paths)

        # chrom sizes and bins
        self.chrom_sizes_path = chrom_sizes_path
        self.chrom_sizes = read_chromsizes(self.chrom_sizes_path)
        if bin_size is None:
            self.bin_size = self._get_bin_size_from_cool_ds()
        else:
            self.bin_size = bin_size
        self.bins_df = binnify(self.chrom_sizes, binsize=self.bin_size)
        self.chrom_offsets = _get_chrom_offsets(self.bins_df)

        # sample weights
        self.sample_dim = sample_dim
        self.sample_weights = sample_weights
        if sample_weights is not None:
            self.sample_weights.index.name = self.sample_dim

        self._chrom_ds_cache = {}

        return

    def _get_bin_size_from_cool_ds(self):
        _path = self.cool_ds_paths[0]
        _chrom = self.chrom_sizes.index[0]
        try:
            with open(f"{_path}/{_chrom}/.zattrs") as f:
                attrs = json.load(f)
            bin_size = attrs["cooler_bin_size"]
            return bin_size
        except BaseException as e:
            print("Cannot infer bin size from cool ds")
            raise e

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
                decode_cf=False,
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

    def get_cooler(
        self,
        samples,
        value_type,
        output_prefix,
        da_name="real",
        dtype="float32",
        scale_factor=1,
        zoomify=True,
        zoomify_cpu=1,
        cooler_kwargs=None,
    ):
        """
        Get cooler from cool ds

        Parameters
        ----------
        samples :
            Samples to fetch, if multiple samples are provided, the matrix is summed with sample weights
        value_type :
            Value type to fetch
        output_prefix :
            Output prefix, suffix will be ".cool" if zoomify is False, or ".mcool" if zoomify is True
        da_name :
            Data array name in xarray dataset
        dtype :
            Data type of the returned matrix
        scale_factor :
            Scale factor to apply to matrix
        zoomify :
            Zoomify the matrix
        zoomify_cpu :
            Number of CPUs to use for zoomify
        cooler_kwargs :
            Additional arguments to pass to create_cooler
        """
        with dask.config.set(scheduler="sync"):
            import subprocess

            from cooler import create_cooler
            from scipy.sparse import coo_matrix

            def _chrom_iterator(
                _samples,
                _value_type,
                _chrom_offset,
                _da_name,
                add_trans=False,
            ):
                """Iterate through the raw matrices and chromosomes of cells."""
                chrom_sizes = self.chrom_sizes

                def _iter_1d(_chrom1, _chrom2):
                    print(f"Saving {_chrom1} x {_chrom2 if _chrom2 is not None else _chrom1}...")
                    chrom_ds = self.fetch(chrom=_chrom1, chrom2=_chrom2)

                    # get chrom 2D np.array
                    matrix = chrom_ds.matrix(
                        samples=_samples,
                        value_type=_value_type,
                        scale_factor=scale_factor,
                        fill_lower_triangle=False,
                        log1p=False,
                        da_name=_da_name,
                        rotate=False,
                        rotate_cval=np.NaN,
                        rotate_height_bp=5000000,
                        dtype=dtype,
                    )
                    if len(matrix.shape) > 2:
                        matrix = matrix.squeeze()

                    # to coo then to pixel
                    matrix = coo_matrix(matrix)
                    _pixel_df = pd.DataFrame({"bin1_id": matrix.row, "bin2_id": matrix.col, "count": matrix.data})

                    # add chrom offset
                    if _chrom2 is None:
                        # both row and col are chrom1
                        _pixel_df.iloc[:, :2] += _chrom_offset[_chrom1]
                    else:
                        # row is chrom1, add chrom1 offset
                        _pixel_df.iloc[:, 0] += _chrom_offset[_chrom1]
                        # col is chrom2, add chrom2 offset
                        _pixel_df.iloc[:, 1] += _chrom_offset[_chrom2]
                    return _pixel_df

                if add_trans:
                    raise NotImplementedError
                else:
                    for chrom in chrom_sizes.keys():
                        pixel_df = _iter_1d(chrom, None)
                        yield pixel_df

            bins_df = binnify(self.chrom_sizes, binsize=self.bin_size)
            chrom_offset = _get_chrom_offsets(bins_df)

            create_cooler(
                cool_uri=f"{output_prefix}.cool",
                bins=bins_df,
                pixels=_chrom_iterator(
                    _samples=samples,
                    _value_type=value_type,
                    _chrom_offset=chrom_offset,
                    _da_name=da_name,
                    add_trans=False,
                ),
                dtypes={"count": dtype},
                ordered=True,
                **(cooler_kwargs or {}),
            )

            if zoomify:
                subprocess.run(["cooler", "zoomify", f"{output_prefix}.cool", "-p", str(zoomify_cpu)], check=True)
                # delete the original cooler file
                subprocess.run(["rm", f"{output_prefix}.cool"], check=True)

        return


class CoolDSChrom(xr.Dataset):
    __slots__ = ()

    def __init__(self, dataset, sample_dim, sample_weights, bin_size):
        super().__init__(data_vars=dataset.data_vars, coords=dataset.coords, attrs=dataset.attrs)
        self.attrs["sample_dim"] = sample_dim

        if bin_size is not None:
            self.attrs["bin_size"] = bin_size

        if sample_weights is not None:
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
        rotate: bool
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
        with dask.config.set(scheduler="synchronous"):
            sel_dict = {}
            if samples is not None:
                sel_dict[self.sample_dim] = samples
            if value_type is not None:
                if da_name not in self.data_vars:
                    value_type_dim_name = "value_type"
                else:
                    value_type_dim_name = f"{da_name}_value_type"
                sel_dict[value_type_dim_name] = value_type
            if len(sel_dict) > 0:
                sample_da = self.sel(sel_dict)[da_name]
            else:
                sample_da = self[da_name]

            if samples is not None and not isinstance(samples, str):
                # sum sample_dim if multiple samples selected
                use_weights = self.sample_weights.sel({self.sample_dim: samples})
                sample_da = (sample_da * use_weights).sum(dim=self.sample_dim) / use_weights.sum() * scale_factor

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
