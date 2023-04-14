import pathlib
import subprocess

import numpy as np
import pandas as pd
import pyBigWig
import xarray as xr
from numba import njit

from ALLCools.utilities import parse_mc_pattern


def _mvalue(mc, cov, alpha=0.1):
    m = np.log2((mc + alpha) / (cov - mc + alpha))
    return m


def _frac(mc, cov, alpha=0.001):
    frac = mc / (cov + alpha)
    return frac


@njit
def _regions_to_pos(regions: np.ndarray) -> np.ndarray:
    """Convert regions to positions."""
    delta = regions[:, 1] - regions[:, 0]
    pos_sel = np.zeros(delta.sum(), dtype=np.uint32)
    cur_pos = 0
    for start, end in regions:
        length = end - start
        pos_sel[cur_pos : cur_pos + length] = np.arange(start, end, dtype=np.uint32)
        cur_pos += length
    return pos_sel


def _chunk_pos_to_bed_df(chrom, chunk_pos):
    """Convert chunk positions to bed dataframe."""
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

        if mc_pattern is None:
            # get all mc types
            judge = np.ones_like(self.mc_type, dtype=bool)
        else:
            # get mc types matching the pattern
            judge = self.mc_type.isin(parse_mc_pattern(mc_pattern))
        # value can be -1, 0, 1, only 0 is False, -1 and 1 are True
        _bool = self.sel(mc_type=judge).sum(dim="mc_type").values.astype(bool)
        return _bool

    def get_mc_pos(self, mc_pattern, offset=None):
        """Get the positions of mc types matching the pattern."""
        mc_pattern = self._validate_mc_pattern(mc_pattern)

        _bool = self.get_mc_pos_bool(mc_pattern)

        _pos = self.get_index("pos")[_bool].copy()
        if offset is not None:
            _pos += offset
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
        return

    @property
    def continuous(self):
        return self.attrs.get("continuous", False)

    @continuous.setter
    def continuous(self, value):
        assert isinstance(value, bool), "continuous must be a boolean."
        self.attrs["continuous"] = value
        if not value:
            # remove offset when continuous is False
            self.attrs.pop("offset", None)

    @property
    def offset(self):
        value = self.attrs.get("offset", 0)
        return int(value)

    @offset.setter
    def offset(self, value):
        assert isinstance(value, int), "offset must be an integer."
        self.attrs["offset"] = value

    def fetch(self, start=None, end=None):
        """
        Select by single region.

        Coordinates are 0-based and half-open, like the BED format.

        Parameters
        ----------
        start :
            The start position of region.
        end :
            The end position of region.

        Returns
        -------
        BaseDSChrom
        """
        if start is None:
            start = 0
        if end is None:
            end = self.chrom_size

        assert start < end, "start must be less than end."

        if self.continuous:
            obj = self.sel(pos=slice(start - self.offset, end - self.offset))
            # copy attrs to avoid changing the original attrs of self
            obj.attrs = obj.attrs.copy()
            obj.attrs["offset"] = start
        else:
            obj = self.fetch_regions([(start, end)])
            # attrs will be copied in fetch_regions
        return obj

    def fetch_regions(self, regions):
        """
        Select by multiple regions.

        Parameters
        ----------
        regions :
            The regions to select. Iterable of (start, end) tuples.
            Coordinates are 0-based and half-open, like the BED format.

        Returns
        -------
        BaseDSChrom
        """
        regions = np.array(regions).astype("uint32")

        assert regions.shape[1] == 2, "Regions can not be converted to (N, 2) array."

        pos_sel = _regions_to_pos(regions)
        obj = self.fetch_positions(pos_sel)
        # attrs will be copied in fetch_positions

        assert obj.continuous is False
        assert obj.offset == 0
        return obj

    def _fetch_positions(self, positions):
        """Fetch positions only, no prefetch."""
        if positions.size == 0:
            obj = self.isel(pos=[])
        else:
            positions = np.sort(positions).astype("uint32")

            # always apply offset to positions
            # if continuous is False, offset will be 0, so no effect
            obj = self.sel(pos=positions - self.offset)

        # copy attrs to avoid changing the original attrs of self
        obj.attrs = obj.attrs.copy()
        obj.continuous = False
        obj.coords["pos"] = positions
        return obj

    def fetch_positions(self, positions):
        """
        Select the positions to create a discontinuous BaseDSChrom.

        Coordinates are 0-based.

        Parameters
        ----------
        positions :
            using genome position to select the positions

        Returns
        -------
        BaseDSChrom
        """
        if self.continuous and len(positions) != 0:
            # if obj is continuous, pre-fetch the min max pos to speed up
            pos_min = min(positions)
            pos_max = max(positions) + 1
            obj = self.fetch(pos_min, pos_max)
            return obj._fetch_positions(positions)
        else:
            return self._fetch_positions(positions)

    @staticmethod
    def _xarray_open(path, obs_dim):
        multi = False
        if isinstance(path, (str, pathlib.Path)):
            if "*" in str(path):
                multi = True
        else:
            if len(path) > 1:
                multi = True
            else:
                path = path[0]

        if multi:
            ds = xr.open_mfdataset(path, concat_dim=obs_dim, combine="nested", engine="zarr", decode_cf=False)
        else:
            ds = xr.open_zarr(path, decode_cf=False)
        return ds

    @classmethod
    def open(cls, path, start=None, end=None, codebook_path=None, obs_dim="sample_id"):
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
        obs_dim
            The dimension name of the observation dimension.

        Returns
        -------
        BaseDSChrom
        """
        _zarr_obj = cls._xarray_open(path, obs_dim=obs_dim)

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
        _obj.continuous = True

        if start is not None or end is not None:
            _obj = _obj.fetch(start, end)
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
        cb = Codebook(self["codebook"])
        cb.attrs["c_pos"] = self.attrs["c_pos"]
        cb.attrs["context_size"] = self.attrs["context_size"]
        return cb

    @property
    def cb(self) -> Codebook:
        """Alias for codebook."""
        return self.codebook

    def select_mc_type(self, pattern):
        """
        Select the mc_type by pattern.

        Parameters
        ----------
        pattern :
            The pattern to select the mc_type.

        Returns
        -------
        BaseDSChrom
        """
        pattern_pos = self.codebook.get_mc_pos(pattern, offset=self.offset)

        ds = self.fetch_positions(positions=pattern_pos)
        assert ds.continuous is False
        return ds

    @property
    def pos_index(self):
        """The position index."""
        return self.get_index("pos")

    def get_region_ds(
        self,
        mc_type,
        bin_size=None,
        regions=None,
        region_name=None,
        region_chunks=10000,
        region_start=None,
        region_end=None,
    ):
        """
        Get the region dataset.

        Parameters
        ----------
        mc_type
            The mc_type to be selected.
        bin_size
            The bin size to aggregate BaseDS to fix-sized regions.
        regions
            The regions dataframe containing tow (start, end) or three columns (chrom, start, end).
            The index will be used as the region names.
        region_name
            The dimension name of the regions.
        region_chunks
            The chunk size of the region dim in result dataset.
        region_start
            The start position of the region to be selected.
            Relevant only when bin_size is provided and regions is None.
        region_end
            The end position of the region to be selected.
            Relevant only when bin_size is provided and regions is None.

        Returns
        -------
        BaseDSChrom
        """
        if bin_size is None and regions is None:
            raise ValueError("One of bin_size or regions must be specified.")

        if bin_size is not None:
            assert bin_size > 1, "bin_size must be greater than 1."
            all_idx = self.get_index("pos")
            region_start = all_idx.min() + self.offset if region_start is None else region_start
            region_end = all_idx.max() + self.offset if region_end is None else region_end

        if regions is not None:
            if regions.shape[1] == 3:
                _chrom_sel = regions.iloc[:, 0] == self.chrom
                regions = regions.loc[_chrom_sel,].iloc[:, 1:].copy()
            assert regions.shape[1] == 2
            assert regions.shape[0] > 0, "regions must have at least one row."

        # get positions
        if mc_type is not None:
            pos_idx = self.cb.get_mc_pos(mc_type, offset=self.offset)
            if regions is not None:
                region_pos_idx = _regions_to_pos(regions=np.array(regions))
                pos_idx = pos_idx.intersection(region_pos_idx).astype(int)
            base_ds = self.fetch_positions(positions=pos_idx)
        else:
            # use all positions in the BaseDS (mc_type selection may have been done)
            base_ds = self
            pos_idx = base_ds.get_index("pos")
            if regions is not None:
                region_pos_idx = _regions_to_pos(regions=np.array(regions))
                pos_idx = pos_idx.intersection(region_pos_idx).astype(int)
                base_ds = self.fetch_positions(positions=pos_idx)
                # update pos_idx
                pos_idx = base_ds.get_index("pos")

        # prepare regions
        if regions is not None:
            bins = regions.iloc[:, 0].tolist() + [regions.iloc[-1, 1]]
            labels = regions.index
            region_name = regions.index.name
        else:
            bins = []
            for i in range(0, self.chrom_size, bin_size):
                if i < region_start or i > region_end:
                    continue
                bins.append(i)
            if bins[-1] < region_end:
                bins.append(region_end)
            if bins[0] > region_start:
                bins.insert(0, region_start)

            labels = []
            for start in bins[:-1]:
                labels.append(start)

        if mc_type is not None and len(pos_idx) == 0:
            # border case, no CpX selected by mc_type
            # create an empty region_ds
            region_ds = base_ds["data"].rename({"pos": "pos_bins"})
            region_ds = region_ds.reindex({"pos_bins": labels}, fill_value=0)
        else:
            region_ds = (
                base_ds["data"]
                .groupby_bins(
                    group="pos",
                    bins=bins,
                    # important, include the start position but exclude the end position;
                    # the same as BED format
                    right=False,
                    labels=labels,
                    precision=3,
                    include_lowest=True,
                    squeeze=True,
                    restore_coord_dims=False,
                )
                .sum(dim="pos")
            )
        if region_chunks is not None:
            region_ds = region_ds.chunk({"pos_bins": region_chunks})
        if region_name is not None:
            region_ds = region_ds.rename({"pos_bins": region_name})

        region_ds = xr.Dataset({"data": region_ds}).fillna(0)
        return region_ds

    def call_dms(
        self,
        groups,
        output_path=None,
        mcg_pattern="CGN",
        n_permute=3000,
        alpha=0.01,
        max_row_count=50,
        max_total_count=3000,
        filter_sig=True,
        merge_strand=True,
        estimate_p=True,
        cpu=1,
        **output_kwargs,
    ):
        """
        Call DMS for a genomic region.

        Parameters
        ----------
        groups :
            Grouping information for the samples.
            If None, perform DMS test on all samples in the BaseDS.
            If provided, first group the samples by the group information, then perform DMS test on each group.
            Samples not occur in the group information will be ignored.
        output_path :
            Path to the output DMS dataset.
            If provided, the result will be saved to disk.
            If not, the result will be returned.
        mcg_pattern :
            Pattern of the methylated cytosine, default is "CGN".
            If None, use all positions in the BaseDS.
        n_permute :
            Number of permutation to perform.
        alpha :
            Minimum p-value/q-value to consider a site as significant.
        max_row_count :
            Maximum number of base counts for each row (sample) in the DMS input count table.
        max_total_count :
            Maximum total number of base counts in the DMS input count table.
        estimate_p :
            Whether to estimate p-value by approximate the null distribution of S as normal distribution.
            The resolution of the estimated p-value is much higher than the exact p-value,
            which is necessary for multiple testing correction.
            FDR corrected q-value is also estimated if this option is enabled.
        filter_sig :
            Whether to filter out the non-significant sites in output DMS dataset.
        merge_strand :
            Whether to merge the base counts of CpG sites next to each other.
        cpu :
            Number of CPU to use.
        output_kwargs :
            Keyword arguments for the output DMS dataset, pass to xarray.Dataset.to_zarr.

        Returns
        -------
        xarray.Dataset if output_path is None, otherwise None.
        """
        from ..dmr.call_dms_baseds import call_dms_worker

        # TODO validate if the BaseDS has the required data for calling DMS

        dms_ds = call_dms_worker(
            groups=groups,
            base_ds=self,
            mcg_pattern=mcg_pattern,
            n_permute=n_permute,
            alpha=alpha,
            max_row_count=max_row_count,
            max_total_count=max_total_count,
            estimate_p=estimate_p,
            cpu=cpu,
            chrom=self.chrom,
            filter_sig=filter_sig,
            merge_strand=merge_strand,
            output_path=output_path,
            **output_kwargs,
        )
        return dms_ds

    def get_sample_mc_values(self, da_name, sample_dim, sample, value):
        da = self[da_name]
        mc = da.sel({sample_dim: sample, "count_type": "mc"})
        cov = da.sel({sample_dim: sample, "count_type": "cov"})
        if value == "frac":
            values = _frac(mc=mc, cov=cov).values
        elif value == "mvalue":
            values = _mvalue(mc=mc, cov=cov).values
        else:
            raise ValueError(f"Unknown value: {value}")
        return values

    def dump_sample_bigwig(self, bigwig_path, da_name, sample_dim, sample, value="frac"):
        pos = self.get_index("pos")
        values = self.get_sample_mc_values(da_name=da_name, sample_dim=sample_dim, sample=sample, value=value)

        # save to bigwig
        bigwig_path = pathlib.Path(bigwig_path).absolute().resolve()
        # noinspection PyArgumentList
        with pyBigWig.open(str(bigwig_path), "w") as bw:
            bw.addHeader([(self.chrom, self.chrom_size)])
            bw.addEntries(self.chrom, pos.values, values=values, span=1)
        return bigwig_path

    def dump_multi_sample_bigwig(self, bigwig_dir, sample_dim, da_name, samples=None, value="frac"):
        sample_idx = self.get_index(sample_dim)
        bigwig_dir = pathlib.Path(bigwig_dir).resolve().absolute()
        bigwig_dir.mkdir(parents=True, exist_ok=True)

        bigwig_paths = []
        if samples is not None:
            sample_idx = sample_idx.intersection(samples)

        for sample in sample_idx:
            bigwig_path = bigwig_dir / f"{sample_dim}-{sample}-{value}.bw"

            self.dump_sample_bigwig(
                bigwig_path=bigwig_path, da_name=da_name, sample_dim=sample_dim, sample=sample, value=value
            )
            bigwig_paths.append(bigwig_path)
        return bigwig_paths

    def dump_base_position_bed(self, bed_path, use_pos=None, bgzip_tabix=True, merge_cpg_strand=False):
        bed_path = str(bed_path)
        if bed_path.endswith(".gz"):
            bed_path = bed_path[:-3]
        bed_path = pathlib.Path(bed_path).resolve().absolute()

        pos = self.get_index("pos")
        if use_pos is not None:
            if isinstance(use_pos, str):
                use_pos = self[use_pos].values
                pos = pos[use_pos]
            else:
                pos = pos.intersection(use_pos)

        if merge_cpg_strand:
            # applicable for CpG only
            pos_and_strand = self["codebook"].load().sum(dim="mc_type").to_pandas()
            forward_pos = pos_and_strand[pos_and_strand > 0].index
            pos = pos.intersection(forward_pos)
            pos_delta = 2
        else:
            pos_delta = 1

        bed_df = pd.DataFrame({"chrom": self.chrom, "start": pos, "end": pos + pos_delta})
        # noinspection PyTypeChecker
        bed_df.to_csv(bed_path, sep="\t", header=None, index=None)
        if bgzip_tabix:
            subprocess.run(f"bgzip -f {bed_path} && " f"tabix -p bed -f {bed_path}.gz", shell=True, check=True)
        return bed_path

    def create_jbrowse(
        self,
        sample_dim,
        da_name="data",
        samples=None,
        bigwig_value="frac",
        add_pos_track=False,
        use_pos=None,
        merge_cpg_strand=False,
        jbrowse_path=None,
        config=None,
        fasta_path=None,
        gene_gtf=None,
        transcript_gtf=None,
        verbose=False,
    ):
        from ..jbrowse import JBrowse

        jb = JBrowse(path=jbrowse_path, config=config, verbose=verbose)

        # add genome reference tracks
        if fasta_path is not None or gene_gtf is not None or transcript_gtf is not None:
            jb.create(fasta_path=fasta_path, gene_gtf=gene_gtf, transcript_gtf=transcript_gtf)

        # add BaseDS tracks
        bigwig_paths = self.dump_multi_sample_bigwig(
            bigwig_dir="/tmp/", sample_dim=sample_dim, da_name=da_name, samples=samples, value=bigwig_value
        )
        for bigwig_path in bigwig_paths:
            # bigwig is already in the jbrowse dir
            jb.add_track(track_path=bigwig_path, name=None, load="move")

        # add position regions
        if add_pos_track:
            self.dump_base_position_bed(
                bed_path="/tmp/BaseDS.Positions.bed.gz",
                use_pos=use_pos,
                bgzip_tabix=True,
                merge_cpg_strand=merge_cpg_strand,
            )
            jb.add_track(track_path="/tmp/BaseDS.Positions.bed.gz", name="BaseDS Positions", load="move")
        return jb


class BaseDS:
    def __init__(self, paths, chrom_sizes_path, codebook_path=None):
        """
        A wrapper for one or multiple BaseDS datasets.

        Parameters
        ----------
        paths :
            Path to the BaseDS datasets.
        chrom_sizes_path :
            Path to the chromosome sizes file.
        codebook_path :
            Path to the codebook file.
        """
        self.paths = self._parse_paths(paths)
        self.codebook_path = codebook_path
        self.chrom_sizes_path = chrom_sizes_path
        self.chrom_sizes = pd.read_csv(chrom_sizes_path, sep="\t", header=None, index_col=0).squeeze()
        self.__base_ds_cache = {}

    @staticmethod
    def _parse_paths(paths):
        import glob

        _paths = []
        if isinstance(paths, str):
            if "*" in paths:
                _paths += list(glob.glob(paths))
            else:
                _paths.append(paths)
        elif isinstance(paths, pathlib.Path):
            _paths.append(str(paths))
        else:
            _paths += list(paths)
        return _paths

    def _get_chrom_paths(self, chrom):
        return [f"{p}/{chrom}" for p in self.paths]

    def _get_chrom_ds(self, chrom):
        if chrom not in self.__base_ds_cache:
            self.__base_ds_cache[chrom] = BaseDSChrom.open(
                path=self._get_chrom_paths(chrom),
                codebook_path=f"{self.codebook_path}/{chrom}",
            )
        _chrom_ds: BaseDSChrom = self.__base_ds_cache[chrom]
        return _chrom_ds

    def fetch(self, chrom, start=None, end=None):
        """
        Fetch a BaseDS for a genomic region.

        Parameters
        ----------
        chrom :
            Chromosome name.
        start :
            Select genomic region by start position.
        end :
            Select genomic region by end position.

        Returns
        -------
        BaseDSChrom
        """
        _chrom_ds = self._get_chrom_ds(chrom)
        return _chrom_ds.fetch(start=start, end=end)

    def fetch_regions(self, chrom, regions):
        """
        Fetch a BaseDS for a list of genomic regions.

        Parameters
        ----------
        chrom :
            Chromosome name.
        regions :
            The regions to select. Iterable of (start, end) tuples.
            Coordinates are 0-based and half-open, like the BED format.

        Returns
        -------
        BaseDSChrom
        """
        _chrom_ds = self._get_chrom_ds(chrom)
        return _chrom_ds.fetch_regions(regions)

    def fetch_positions(self, chrom, positions):
        """
        Fetch a BaseDS for a list of genomic positions.

        Coordinates are 0-based.

        Parameters
        ----------
        chrom :
            Chromosome name.
        positions :
            using genome position to select the positions

        Returns
        -------
        BaseDSChrom
        """
        _chrom_ds = self._get_chrom_ds(chrom)
        return _chrom_ds.fetch_positions(positions)
