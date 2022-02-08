import pandas as pd
import pybedtools
import xarray as xr
import numpy as np
import dask
import warnings
import joblib
import subprocess
import pathlib
import yaml
import pyBigWig
from pybedtools import BedTool
from concurrent.futures import ProcessPoolExecutor, as_completed
from .utilities import determine_engine, obj_to_str, write_ordered_chunks, update_dataset_config
import os
from ALLCools.utilities import parse_chrom_size

os.environ["NUMEXPR_MAX_THREADS"] = "16"


def _bigwig_over_bed(bed: pd.DataFrame, path, value_type="mean", dtype="float32"):
    with pyBigWig.open(path, "r") as bw:

        def _region_stat(row, t=value_type):
            chrom, start, end, *_ = row
            try:
                value = bw.stats(chrom, start, end, type=t)[0]
            except RuntimeError:
                # happens when the region has error or chrom not exist in bw
                # let user decide what happen, here just return nan
                value = np.NaN
            return value

        values = bed.apply(_region_stat, t=value_type, axis=1)
        values = values.astype(dtype)
    return values


def _region_bed_sorted(bed_path, g, bed_sorted):
    chrom_sizes = parse_chrom_size(g)

    bed_df = pd.read_csv(bed_path, sep="\t", index_col=None, header=None)
    # select chroms that exist in g
    bed_df = bed_df.loc[bed_df.iloc[:, 0].isin(chrom_sizes.keys())]
    bed = BedTool.from_dataframe(bed_df)

    if bed_sorted:
        return bed
    else:
        return bed.sort(g=g)


def _bed_intersection(bed: pybedtools.BedTool, path, g, region_index, bed_sorted, fraction=0.2):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        query_bed = _region_bed_sorted(path, g, bed_sorted)
        try:
            df = bed.intersect(
                query_bed, wa=True, f=fraction, g=g, sorted=True
            ).to_dataframe()
            if df.shape[0] == 0:
                regions_idx = pd.Series([])
            else:
                regions_idx = df["name"]
        except pd.errors.EmptyDataError:
            regions_idx = pd.Series([])
    regions = pd.Index(regions_idx.values)
    bool_series = pd.Series(region_index.isin(regions), index=region_index)

    query_bed.delete_temporary_history(ask=False)
    return bool_series


def _annotate_by_bigwigs_worker(
        dataset_path,
        region_dim,
        chrom_size_path,
        track_paths,
        output_path,
        dim,
        slop,
        value_type,
        dtype,
        **kwargs,
):
    len(kwargs)
    # set dask scheduler to allow multiprocessing
    with dask.config.set(scheduler="sync"):
        # Open region ds again inside the worker function
        region_ds = RegionDS.open(
            path=dataset_path, region_dim=region_dim, chrom_size_path=chrom_size_path
        )

        # get dmr region bed and bigwig files
        dmr_bed = region_ds.get_bed(
            with_id=False, bedtools=False, slop=slop, chrom_size_path=chrom_size_path
        )
        # iterate each bigwig
        total_values = {}
        for sample, bigwig_path in track_paths.items():
            values = _bigwig_over_bed(
                bed=dmr_bed, path=bigwig_path, value_type=value_type, dtype=dtype
            )
            total_values[sample] = values
        total_values = pd.DataFrame(total_values)
        total_values.columns.name = dim
        total_values.index.name = region_dim

        ds = xr.Dataset({f"{region_dim}_{dim}_da": total_values})
        ds.to_zarr(output_path, mode="w")
    return output_path


def _annotate_by_beds_worker(
        dataset_path,
        region_dim,
        chrom_size_path,
        slop,
        track_paths,
        dtype,
        dim,
        output_path,
        bed_sorted,
        fraction=0.2,
        **kwargs,
):
    len(kwargs)

    # set dask scheduler to allow multiprocessing
    with dask.config.set(scheduler="sync"):
        # Open region ds again inside the worker function
        region_ds = RegionDS.open(
            path=dataset_path, region_dim=region_dim, chrom_size_path=chrom_size_path
        )

        # get dmr region bed
        dmr_bed = region_ds.get_bed(
            with_id=True, bedtools=True, slop=slop, chrom_size_path=chrom_size_path
        ).sort(g=chrom_size_path)

        total_values = {}
        for sample, bed_path in track_paths.items():
            values = _bed_intersection(
                bed=dmr_bed,
                path=bed_path,
                bed_sorted=bed_sorted,
                g=chrom_size_path,
                region_index=region_ds.get_index(region_ds.region_dim),
                fraction=fraction
            )
            total_values[sample] = values.astype(dtype)
        total_values = pd.DataFrame(total_values)
        total_values.columns.name = dim

        ds = xr.Dataset({f"{region_dim}_{dim}_da": total_values})
        ds.to_zarr(output_path, mode="w")

        dmr_bed.delete_temporary_history(ask=False)
    return output_path


def _fisher_exact(row, alternative="two-sided"):
    from scipy.stats import fisher_exact

    oddsratio, p = fisher_exact(row.values.reshape((2, 2)), alternative=alternative)
    value = pd.Series({"oddsratio": oddsratio, "p": p})
    return value


class RegionDS(xr.Dataset):
    __slots__ = ()

    def __init__(self, dataset, region_dim=None, location=None, chrom_size_path=None):
        super().__init__(
            data_vars=dataset.data_vars, coords=dataset.coords, attrs=dataset.attrs
        )
        self.region_dim = region_dim
        self.location = location
        self.chrom_size_path = chrom_size_path
        return

    @property
    def region_dim(self):
        return self.attrs.get("region_dim")

    @region_dim.setter
    def region_dim(self, region_dim):
        if region_dim is not None:
            if region_dim not in self.dims:
                raise KeyError(
                    f"{region_dim} does not occur in dimension names: {list(self.dims.keys())}"
                )
            self.attrs["region_dim"] = region_dim
        else:
            return

    @property
    def chrom_size_path(self):
        return self.attrs.get("chrom_size_path")

    @chrom_size_path.setter
    def chrom_size_path(self, chrom_size_path):
        if chrom_size_path is not None:
            chrom_size_path = pathlib.Path(chrom_size_path).absolute()
            if not chrom_size_path.exists():
                raise FileNotFoundError(str(chrom_size_path))
            self.attrs["chrom_size_path"] = str(chrom_size_path)
        else:
            return

    @property
    def location(self):
        return self.attrs.get("region_ds_location")

    @location.setter
    def location(self, path):
        if path is not None:
            location = pathlib.Path(path).absolute()
            self.attrs["region_ds_location"] = str(location)
            location.mkdir(exist_ok=True, parents=True)
        else:
            return

    @classmethod
    def from_bed(
            cls, bed, location, chrom_size_path, region_dim="region", sort_bed=True
    ):
        """
        Create empty RegionDS from a bed file.

        Parameters
        ----------
        bed
        location
        region_dim
        chrom_size_path
        sort_bed

        Returns
        -------

        """

        # sort bed based on chrom_size_path
        if isinstance(bed, (str, pathlib.PosixPath)):
            if sort_bed:
                bed = BedTool(bed).sort(g=chrom_size_path).to_dataframe()
            else:
                bed = BedTool(bed)
        else:
            bed = bed

        n_cols = bed.shape[1]
        if n_cols == 3:
            bed.index = bed.index.map(lambda i: f"{region_dim}_{i}")
        elif n_cols == 4:
            bed.set_index(bed.columns[3], inplace=True)
        else:
            raise ValueError(
                "bed file need to be either 3 columns (chrom, start, end) "
                "or 4 columns (chrom, start, end, name)"
            )
        bed.index.name = region_dim
        bed.columns = ["chrom", "start", "end"]

        ds = xr.Dataset({})
        region_dim = bed.index.name
        for k, v in bed.items():
            key = f"{region_dim}_{k}"
            ds.coords[key] = v
            if ds.coords[key].dtype == "object":
                ds.coords[key] = ds.coords[key].astype(str)

        location = pathlib.Path(location).absolute()
        location.mkdir(exist_ok=True, parents=True)
        region_ds = cls(
            ds,
            region_dim=region_dim,
            location=location,
            chrom_size_path=chrom_size_path,
        )
        region_ds.save()
        return region_ds

    @classmethod
    def open(
            cls,
            path,
            region_dim=None,
            use_regions=None,
            split_large_chunks=True,
            chrom_size_path=None,
            select_dir=None,
            engine="zarr",
    ):
        if isinstance(path, (str, pathlib.PosixPath)):
            _path = pathlib.Path(path).absolute()

            if _path.is_dir():
                # check if this is a RegionDS dir that contains multiple datasets
                region_ds_file = _path / ".ALLCools"
                if region_ds_file.exists():
                    with open(region_ds_file) as f:
                        region_ds_config = yaml.load(f, yaml.SafeLoader)
                        if region_dim is None:
                            region_dim = region_ds_config["region_dim"]
                            print(f"Using {region_dim} as region_dim")
                        # only open datasets having the region_dim
                        # other datasets will not be opened
                        # e.g, when region_dim == 'dmr', dms dataset will not be opened.
                        exclude_dir_name = [
                            k
                            for k, v in region_ds_config["ds_region_dim"].items()
                            if v != region_dim
                        ]

                    # add output_dir
                    region_ds_location = pathlib.Path(path).absolute()
                    # add chrom size path if exist
                    if chrom_size_path is None:
                        chrom_size_path = _path / "chrom_sizes.txt"
                        if not chrom_size_path.exists():
                            chrom_size_path = None

                    # read all sub dir as a RegionDS, then merge all the RegionDS together.
                    datasets = []
                    for sub_dir_path in _path.iterdir():
                        if sub_dir_path.is_dir() and sub_dir_path.name[0] != ".":
                            if select_dir is not None:
                                if sub_dir_path.name not in select_dir:
                                    # not in select_dir list, skip
                                    continue
                            if sub_dir_path.name in exclude_dir_name:
                                # the dir do not have region_dim, skip
                                continue
                            if not (sub_dir_path / ".zattrs").exists():
                                # print(f'{sub_dir_path} does not seem to be a zarr storage, skipped')
                                continue
                            try:
                                datasets.append(
                                    cls._open_single_dataset(
                                        path=sub_dir_path,
                                        region_dim=region_dim,
                                        split_large_chunks=split_large_chunks,
                                        engine=engine,
                                    )
                                )
                            except BaseException as e:
                                print(f"An error raised when reading {sub_dir_path}.")
                                raise e
                    region_ds = cls(
                        xr.merge(datasets),
                        region_dim=region_dim,
                        location=region_ds_location,
                        chrom_size_path=chrom_size_path,
                    )
                else:
                    # is dir, but is not RegionDS dir, could be just a zarr dir
                    region_ds = cls._open_single_dataset(
                        path=path,
                        region_dim=region_dim,
                        split_large_chunks=split_large_chunks,
                        chrom_size_path=chrom_size_path,
                        engine=engine,
                    )
            else:
                # dataset stored in other format, such as netcdf4
                region_ds = cls._open_single_dataset(
                    path=path,
                    region_dim=region_dim,
                    split_large_chunks=split_large_chunks,
                    chrom_size_path=chrom_size_path,
                    engine=engine,
                )
        else:
            # could be a list of paths, open it as a single dataset
            path = list(path)
            region_ds = cls._open_single_dataset(
                path=path,
                region_dim=region_dim,
                split_large_chunks=split_large_chunks,
                chrom_size_path=chrom_size_path,
                engine=engine,
            )

        if use_regions is not None:
            region_ds = region_ds.sel({region_dim: use_regions})
        return region_ds

    @classmethod
    def _open_single_dataset(
            cls,
            path,
            region_dim,
            split_large_chunks=True,
            chrom_size_path=None,
            location=None,
            engine=None,
    ):
        """
        Take one or multiple RegionDS file paths and create single RegionDS concatenated on region_dim

        Parameters
        ----------
        path
            Single RegionDS path or RegionDS path pattern with wildcard or RegionDS path list
        region_dim
            Dimension name of regions
        split_large_chunks
            Split large dask array chunks if true
        Returns
        -------
        RegionDS
        """
        # print('opening', path)
        if region_dim is None:
            raise ValueError(
                "Please specify a region_dim name when open a normal xr.Dataset with RegionDS."
            )

        if engine is None:
            engine = determine_engine(path)
        # if engine is None:
        #     print(f'Open RegionDS with netcdf4 engine.')
        # else:
        #     print(f'Open RegionDS with {engine} engine.')
        try:
            if (isinstance(path, str) and "*" not in path) or isinstance(
                    path, pathlib.PosixPath
            ):
                ds = xr.open_dataset(path, engine=engine)
            else:
                with dask.config.set(
                        **{"array.slicing.split_large_chunks": split_large_chunks}
                ):
                    if isinstance(path, str):
                        import glob

                        path = sorted([p for p in glob.glob(path)])
                    ds = xr.open_mfdataset(
                        path,
                        parallel=False,
                        combine="nested",
                        concat_dim=region_dim,
                        engine=engine,
                    )
        except Exception as e:
            print(f"Got error when opening {path}")
            print(f"Engine parameter is {engine}")
            raise e
        ds = cls(
            ds,
            region_dim=region_dim,
            location=location,
            chrom_size_path=chrom_size_path,
        ).squeeze()
        return ds

    def iter_index(self, chunk_size=100000, dim=None):
        if dim is None:
            dim = self.region_dim

        index = self.get_index(dim)
        for chunk_start in range(0, index.size, chunk_size):
            use_index = index[chunk_start: chunk_start + chunk_size]
            yield use_index

    def iter_array(self, chunk_size=100000, dim=None, da=None, load=False):
        if dim is None:
            dim = self.region_dim

        if da is None:
            da = f"{dim}_da"
        _da = self[da]

        assert dim in _da.dims

        for _index in self.iter_index(chunk_size=chunk_size, dim=dim):
            use_da = _da.sel({dim: _index})
            if load:
                use_da.load()
            yield use_da

    def get_fasta(
            self,
            genome_fasta,
            output_path,
            slop=None,
            chrom_size_path=None,
            standardize_length=None,
    ):
        bed = self.get_bed(
            with_id=True,
            bedtools=True,
            slop=slop,
            chrom_size_path=chrom_size_path,
            standardize_length=standardize_length,
        )
        bed.getfasta(fo=output_path, nameOnly=True, fi=genome_fasta)
        return

    def get_bed(
            self,
            with_id=True,
            bedtools=False,
            slop=None,
            chrom_size_path=None,
            standardize_length=None,
    ):
        if chrom_size_path is None:
            chrom_size_path = self.chrom_size_path  # will be none if not exist

        region_dim = self.region_dim

        bed_df = pd.DataFrame(
            {
                "chrom": self.coords[f"{region_dim}_chrom"],
                "start": self.coords[f"{region_dim}_start"],
                "end": self.coords[f"{region_dim}_end"],
            }
        )

        # standardize region length, used in motif enrichment analysis
        if standardize_length is not None:
            # standardize_length is an int number
            region_center = bed_df["start"] + (bed_df["end"] - bed_df["start"]) // 2
            bed_df["start"] = region_center - 1
            bed_df["end"] = region_center
            slop = (
                    standardize_length // 2
            )  # use the bedtools slop to extend the center to standard length

        if with_id:
            bed_df["name"] = self.get_index(region_dim).tolist()

        bed = None
        if slop is not None and slop > 0:
            if chrom_size_path is None:
                raise ValueError("Must provide chrom_size_path when slop is not None.")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bed = BedTool.from_dataframe(bed_df).slop(b=slop, g=chrom_size_path)
            if not bedtools:
                bed_df = bed.to_dataframe()

        if bedtools:
            if bed is None:
                bed = BedTool.from_dataframe(bed_df)
            return bed
        else:
            bed_df.index = self.get_index(self.region_dim)
            return bed_df

    def _chunk_annotation_executor(self, annotation_function, cpu, save=True, **kwargs):
        chrom_size_path = kwargs["chrom_size_path"]
        dim = kwargs["dim"]
        chunk_size = kwargs["chunk_size"]
        # add back when submit jobs, use this to control parallel
        track_paths = kwargs.pop("track_paths")

        region_ds_path = self.location
        if region_ds_path is None:
            raise ValueError(f"Must have an on-disk location to annotate bigwigs.")

        if chrom_size_path is None:
            chrom_size_path = self.chrom_size_path
        region_dim = self.region_dim

        # deal with path
        chunk_dir_path = f"{region_ds_path}/.{region_dim}_{dim}_chunks"
        chunk_dir_path = pathlib.Path(chunk_dir_path)
        if chunk_dir_path.exists():
            subprocess.run(f"rm -rf {chunk_dir_path}", shell=True)
        chunk_dir_path.mkdir(exist_ok=True)
        final_dir_path = f"{region_ds_path}/{region_dim}_{dim}"

        n_features = len(track_paths)
        if chunk_size == "auto":
            chunk_size = max(1, n_features // cpu // 2 + 1)
        print(f"Use chunk size {chunk_size}")

        other_kwargs = {
            "region_dim": region_dim,
            "dataset_path": region_ds_path,
            "chrom_size_path": chrom_size_path,
        }
        kwargs.update(other_kwargs)

        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for i, chunk_start in enumerate(range(0, n_features, chunk_size)):
                output_path = f"{chunk_dir_path}/chunk_{i}.zarr"
                kwargs["output_path"] = output_path
                kwargs["track_paths"] = track_paths[
                                        chunk_start: chunk_start + chunk_size
                                        ]
                future = exe.submit(annotation_function, **kwargs)
                futures[future] = i
                # time.sleep(1)

            chunks_to_write = {}
            for i, future in enumerate(as_completed(futures)):
                chunk_i = futures[future]
                output_path = future.result()
                chunks_to_write[chunk_i] = output_path

        if save:
            write_ordered_chunks(
                chunks_to_write=chunks_to_write,
                final_path=final_dir_path,
                append_dim=dim,
                engine="zarr",
                dtype=kwargs["dtype"],
                coord_dtypes=None,
            )
            update_dataset_config(self.location, add_ds_region_dim={f"{region_dim}_{dim}": region_dim})

            # load the newly generated da only
            _ds = xr.open_zarr(final_dir_path)
        else:
            _ds = xr.open_mfdataset(
                [chunks_to_write[k] for k in sorted(chunks_to_write.keys())],
                concat_dim=dim,
                combine="nested",
                engine="zarr",
            ).load()
        self.update(_ds)

        subprocess.run(f"rm -rf {chunk_dir_path}", shell=True)
        return

    def annotate_by_bigwigs(
            self,
            bigwig_table,
            dim,
            slop=100,
            chrom_size_path=None,
            value_type="mean",
            chunk_size="auto",
            dtype="float32",
            cpu=1,
            save=True,
    ):
        if isinstance(bigwig_table, dict):
            track_paths = pd.Series(bigwig_table)
        elif isinstance(bigwig_table, pd.Series):
            track_paths = bigwig_table
        elif isinstance(bigwig_table, str) and bigwig_table.endswith("csv"):
            track_paths = pd.read_csv(
                bigwig_table, index_col=0, squeeze=True, header=None
            )
        else:
            track_paths = pd.read_csv(
                bigwig_table, sep="\t", index_col=0, squeeze=True, header=None
            )

        kwargs = dict(
            track_paths=track_paths,
            dim=dim,
            slop=slop,
            chrom_size_path=chrom_size_path,
            value_type=value_type,
            chunk_size=chunk_size,
            dtype=dtype,
        )
        self._chunk_annotation_executor(
            _annotate_by_bigwigs_worker, cpu=cpu, save=save, **kwargs
        )
        return

    def annotate_by_beds(
            self,
            bed_table,
            dim,
            slop=100,
            chrom_size_path=None,
            chunk_size="auto",
            dtype="bool",
            bed_sorted=True,
            cpu=1,
            fraction=0.2,
            save=True,
    ):
        bed_tmp = pathlib.Path(
            f"./pybedtools_tmp_{np.random.randint(0, 100000)}"
        ).absolute()
        bed_tmp.mkdir(exist_ok=True)
        default_tmp = pybedtools.helpers.get_tempdir()
        pybedtools.helpers.set_tempdir(str(bed_tmp))

        if isinstance(bed_table, dict):
            track_paths = pd.Series(bed_table)
        elif isinstance(bed_table, pd.Series):
            track_paths = bed_table
        elif isinstance(bed_table, str) and bed_table.endswith("csv"):
            track_paths = pd.read_csv(bed_table, index_col=0, squeeze=True, header=None)
        else:
            track_paths = pd.read_csv(
                bed_table, sep="\t", index_col=0, squeeze=True, header=None
            )

        kwargs = dict(
            track_paths=track_paths,
            dim=dim,
            slop=slop,
            chrom_size_path=chrom_size_path,
            chunk_size=chunk_size,
            dtype=dtype,
            bed_sorted=bed_sorted,
            fraction=fraction
        )
        self._chunk_annotation_executor(
            _annotate_by_beds_worker, cpu=cpu, save=save, **kwargs
        )

        subprocess.run(f"rm -rf {bed_tmp}", shell=True)
        # pybedtools actually changed tempfile.tempdir
        # we must put it back otherwise other package might raise problem
        pybedtools.helpers.set_tempdir(default_tmp)
        return

    def get_feature(self, feature_name, dim=None, da_name=None):
        if dim is None:
            try:
                data = self.coords[feature_name].to_pandas()
            except KeyError:
                raise KeyError(
                    f"{feature_name} does not exist in RegionDS.coords, "
                    f"if it belongs to a dataarray, please specify its dim name."
                )
        else:
            if da_name is None:
                da_name = f"{self.region_dim}_{dim}_da"
            data = self[da_name].sel({dim: feature_name}).to_pandas()
        return data

    def scan_motifs(
            self,
            genome_fasta,
            cpu=1,
            standardize_length=500,
            motif_set_path=None,
            chrom_size_path=None,
            combine_cluster=True,
            fnr_fpr_fold=1000,
            chunk_size=None,
            dim="motif",
    ):
        region_ds_path = self.location
        if region_ds_path is None:
            raise ValueError(f"Must have an on-disk location to annotate bigwigs.")

        if chrom_size_path is None:
            chrom_size_path = self.attrs.get("chrom_size_path")
        region_dim = self.region_dim

        fasta_path = f"{region_ds_path}/regions.fasta"
        self.get_fasta(
            genome_fasta,
            output_path=fasta_path,
            slop=None,
            chrom_size_path=chrom_size_path,
            standardize_length=standardize_length,
        )

        from ..motif import MotifSet

        if motif_set_path is not None:
            motif_set: MotifSet = joblib.load(motif_set_path)
        else:
            # default motif database from three sources
            from ..motif import get_default_motif_set

            motif_set = get_default_motif_set()

        if fnr_fpr_fold != 1000:
            motif_set.calculate_threshold(
                cpu=cpu, method="balance", threshold_value=fnr_fpr_fold
            )

        motif_ds = motif_set.scan_motifs(
            fasta_path=fasta_path,
            output_dir=region_ds_path,
            cpu=cpu,
            region_dim=region_dim,
            combine_cluster=combine_cluster,
            motif_dim=dim,
            chunk_size=chunk_size,
        )
        new_dataset_dim = {f'{region_dim}_{dim}': region_dim}
        if combine_cluster:
            new_dataset_dim[f'{region_dim}_{dim}-cluster'] = region_dim
        update_dataset_config(self.location, add_ds_region_dim=new_dataset_dim)
        self.update(motif_ds)
        return

    def get_hypo_hyper_index(
            self,
            a,
            region_dim=None,
            region_state_da=None,
            sample_dim="sample",
            use_collapsed=True,
    ):
        if region_state_da is None:
            if region_dim is None:
                region_dim = self.region_dim
            if use_collapsed:
                region_state_da = f"{region_dim}_state_collapsed"
                if not sample_dim.endswith("_collapsed"):
                    sample_dim = f"{sample_dim}_collapsed"
                if region_state_da not in self.data_vars:
                    print(
                        f'use_collapsed=True but the "{region_dim}_state_collapsed" '
                        f"dataarray do not exist, provide region_state_da name explicitly "
                        f"or set use_collapsed=False"
                    )
                    region_state_da = f"{region_dim}_state"
            else:
                region_state_da = f"{region_dim}_state"
        a_states = self[region_state_da].sel({sample_dim: a}).to_pandas()

        hypo_dmr = a_states == -1
        hypo_dmr = hypo_dmr[hypo_dmr].index

        hyper_dmr = a_states == 1
        hyper_dmr = hyper_dmr[hyper_dmr].index

        return hypo_dmr, hyper_dmr

    def get_pairwise_differential_index(
            self,
            a,
            b,
            dmr_type="hypo",
            region_dim=None,
            region_state_da=None,
            sample_dim="sample",
            use_collapsed=True,
    ):
        a_hypo, a_hyper = self.get_hypo_hyper_index(
            a,
            region_dim=region_dim,
            region_state_da=region_state_da,
            sample_dim=sample_dim,
            use_collapsed=use_collapsed,
        )
        b_hypo, b_hyper = self.get_hypo_hyper_index(
            b,
            region_dim=region_dim,
            region_state_da=region_state_da,
            sample_dim=sample_dim,
            use_collapsed=use_collapsed,
        )

        if dmr_type.lower() == "hypo" or dmr_type == -1:
            # a hypo, b not hypo
            a_not_b = a_hypo[~a_hypo.isin(b_hypo)]
            a_and_b = a_hypo[a_hypo.isin(b_hypo)]
            b_not_a = b_hypo[~b_hypo.isin(a_hypo)]
        elif dmr_type.lower() == "hyper" or dmr_type == 1:
            # a hypo, b not hypo
            a_not_b = a_hyper[~a_hyper.isin(b_hyper)]
            a_and_b = a_hyper[a_hyper.isin(a_hyper)]
            b_not_a = b_hyper[~b_hyper.isin(a_hyper)]
        else:
            raise ValueError(f"Unknown dmr_type {dmr_type}.")

        return a_not_b, a_and_b, b_not_a

    def motif_enrichment(
            self,
            true_regions,
            background_regions,
            region_dim=None,
            motif_dim="motif-cluster",
            motif_da=None,
            alternative="two-sided",
    ):
        if region_dim is None:
            region_dim = self.region_dim
        if motif_da is None:
            motif_da = f"{region_dim}_{motif_dim}_da"
        true_motif = self[motif_da].sel(
            {"motif_value": "n_motifs", region_dim: true_regions}
        )
        true_motif = (true_motif > 0).sum(dim=region_dim).to_pandas()
        true_no_motif = true_regions.size - true_motif

        bkg_motif = self[motif_da].sel(
            {"motif_value": "n_motifs", region_dim: background_regions}
        )
        bkg_motif = (bkg_motif > 0).sum(dim=region_dim).to_pandas()
        bkg_no_motif = background_regions.size - bkg_motif

        contingency_tables = pd.DataFrame(
            {"tp": true_motif, "tn": true_no_motif, "fp": bkg_motif, "fn": bkg_no_motif}
        )

        # add pseudo count
        contingency_tables += 1

        test_results = contingency_tables.apply(
            _fisher_exact, axis=1, alternative=alternative
        )
        from statsmodels.stats.multitest import multipletests

        _, q, *_ = multipletests(test_results["p"], method="fdr_bh")
        test_results["q"] = q

        test_results["log2OR"] = np.log2(test_results["oddsratio"])
        test_results["-lgq"] = -np.log10(q)

        return test_results

    def sample_dmr_motif_enrichment(
            self,
            sample,
            region_dim=None,
            sample_dim="sample",
            motif_dim="motif-cluster",
            region_state_da=None,
            motif_da=None,
            alternative="two-sided",
            use_collapsed=True,
    ):
        hypo_region, hyper_region = self.get_hypo_hyper_index(
            sample,
            region_dim=region_dim,
            region_state_da=region_state_da,
            sample_dim=sample_dim,
            use_collapsed=use_collapsed,
        )
        test_results = self.motif_enrichment(
            hypo_region,
            hyper_region,
            region_dim=region_dim,
            motif_dim=motif_dim,
            motif_da=motif_da,
            alternative=alternative,
        )
        return test_results

    def pairwise_dmr_motif_enrichment(
            self,
            a,
            b,
            dmr_type="hypo",
            region_dim=None,
            region_state_da=None,
            sample_dim="sample",
            motif_dim="motif-cluster",
            motif_da=None,
            alternative="two-sided",
    ):
        a_not_b, _, b_not_a = self.get_pairwise_differential_index(
            a,
            b,
            dmr_type=dmr_type,
            region_dim=region_dim,
            region_state_da=region_state_da,
            sample_dim=sample_dim,
        )
        test_results = self.motif_enrichment(
            a_not_b,
            b_not_a,
            region_dim=region_dim,
            motif_dim=motif_dim,
            motif_da=motif_da,
            alternative=alternative,
        )
        return test_results

    def object_coords_to_string(self, dtypes=None):
        obj_to_str(self, dtypes)
        return

    def save(self, da_name=None, output_path=None, mode="w", change_region_dim=True):
        if self.location is None:
            raise ValueError(f"RegionDS.location is None when trying to save.")
        pathlib.Path(self.location).mkdir(parents=True, exist_ok=True)

        if output_path is None:
            output_path = f"{self.location}/{self.region_dim}"
            update_dataset_config(self.location, add_ds_region_dim={self.region_dim: self.region_dim},
                                  change_region_dim=self.region_dim if change_region_dim else None)

        # turn object coords to fix length string dtype before saving to zarr
        if da_name is None:
            ds_to_save = self
        else:
            ds_to_save = self[[da_name]]
        ds_to_save.object_coords_to_string()
        ds_to_save.to_zarr(output_path, mode=mode)
        return

    def get_coords(self, name):
        return self.coords[name].to_pandas()
